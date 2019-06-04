#include <algorithm>
#include <cmath>
#include <cassert>
#include <iterator> // for distance()
#include <map>
#include <queue>
#include "Utility.h"
#include "Complex.h"
#include "Mesh.h"
#include "SparseMatrix.h"
#include "SectionIntegrals.h"

namespace DDG{

  bool Contains( std::map<FaceCIter, HalfEdgeCIter> &fparent,
		 FaceCIter fi, FaceCIter fj ){
    std::map<FaceCIter, HalfEdgeCIter>::const_iterator i = fparent.find(fi), j = fparent.find(fj);
    return
      ( ( i != fparent.end() && i->second->face == fj ) ||
	( j != fparent.end() && j->second->face == fi ) );
  }

  bool Contains( std::map<VertexCIter, HalfEdgeCIter> &vparent,
		 VertexCIter vi, VertexCIter vj ){
    std::map<VertexCIter, HalfEdgeCIter>::const_iterator i = vparent.find(vi), j = vparent.find(vj);
    return
      ( i != vparent.end() && i->second->vertex == vj ) ||
      ( j != vparent.end() && j->second->vertex == vi );
  }

  // build spanning tree of all triangles; look at co-tree if present
  // to make insertion decisions; if there is a hole someplace add a
  // single boundary edge of that hole to the DualTree
  bool DualTree( std::map<FaceCIter, HalfEdgeCIter> &tree,
		 FaceCIter root,
		 std::map<VertexCIter, HalfEdgeCIter> &cotree ){

    tree.insert( std::pair<FaceCIter, HalfEdgeCIter>( root, root->he ) );

    bool firstBoundary = false; // haven't seen one yet
    // BFS
    std::queue<FaceCIter> Q; Q.push( root );
    while( !Q.empty() ){
      const FaceCIter fi = Q.front(); Q.pop();
      assert( !fi->isBoundary() ); // boundary faces should never be in Q
      // visit all neighbors and make them children
      HalfEdgeCIter he = fi->he;
      do{
	// boundary faces are holes and they don't exist for purposes
	// of this algorithm; note that skipping them here means we
	// are not inserting any dual edges which cross the boundary
	// (of course not! we are building a spanning tree of the
	// triangles!) HOWEVER we make one exception for the first
	// face we encounter
	const FaceCIter fj = he->flip->face;
	if( !fj->isBoundary() ){
	  const VertexCIter vi = he->vertex, vj = he->flip->vertex;
	  if( tree.find(fj) == tree.end() && !Contains( cotree, vi, vj ) ){
	    // not yet in there so link it; the guy on the other side
	    // (first) has the current face (second) to point to
	    tree.insert( std::pair<FaceCIter, HalfEdgeCIter>( fj, he ) );
	    Q.push( fj );
	  }
	}else if( !firstBoundary ){
	  firstBoundary = true;
	  const VertexCIter vi = he->vertex, vj = he->flip->vertex;
	  if( tree.find(fj) == tree.end() && !Contains( cotree, vi, vj ) ){
	    tree.insert( std::pair<FaceCIter, HalfEdgeCIter>( fj, he ) );
	  }
	}
      }while( ( he = he->next ) != fi->he );
    }
    return firstBoundary;
  }

  // here we need the triangle spanning tree to identify primal edges
  // which are not crossed by an edge in the dual tree
  void PrimalTree( std::map<VertexCIter, HalfEdgeCIter> &tree,
		   VertexCIter root,
		   std::map<FaceCIter, HalfEdgeCIter> &cotree ){

    tree.insert( std::pair<VertexCIter, HalfEdgeCIter>( root, root->he ) );

    std::queue<VertexCIter> Q; Q.push( root );
    while( !Q.empty() ){
      const VertexCIter vi = Q.front(); Q.pop();
      HalfEdgeCIter he = vi->he;
      do{
	const VertexCIter vj = he->flip->vertex;
	const FaceCIter fi = he->face, fj = he->flip->face;
	// primal edge not yet linked
	if( tree.find(vj) == tree.end() && !Contains( cotree, fi, fj ) ){
	  // add it
	  tree.insert( std::pair<VertexCIter, HalfEdgeCIter>( vj, he ) );
	  Q.push( vj );
	}
      }while( ( he = he->flip->next ) != vi->he );
    }
  }

  void TraceToRoot( std::vector<HalfEdgeCIter> &c, VertexCIter v,
		    std::map<VertexCIter, HalfEdgeCIter> &vparent ){
    // use this for lookup of parent
    std::map<VertexCIter, HalfEdgeCIter>::const_iterator vp;
    // get parent of v; v better be in the tree; if we are not at the root
    while( ( vp = vparent.find(v), assert( vp != vparent.end() ), vp->second->vertex != v ) ){
      // enter into chain leading to root
      c.push_back( vp->second );
      // and move up the tree
      v = vp->second->vertex;
    }
  }

  bool Mesh::FindGenerators( void ){

    // tree and co-tree are represented as parent pointers; given one
    // type a HalfEdge is returned (if any) which gives the
    // corresponding parent
    std::map<VertexCIter, HalfEdgeCIter> vtree;
    std::map<FaceCIter, HalfEdgeCIter> ftree;

    // now build the tree/co-tree; tree is the DUAL spanning tree (it
    // will also have one boundary edge if one is ever encountered)
    FaceCIter froot = faces.begin(); while( froot->isBoundary() ) froot++;
    const bool boundary = DualTree( ftree, froot, vtree );

    // now the cotree for the vertices. Nothing special here
    PrimalTree( vtree, vertices.begin(), ftree );

    // ok; tree (dual) and co-tree (primal) are built; now look for
    // edges which are in neither
    for( EdgeIter ei = edges.begin(); ei != edges.end(); ei++ ){
      if( !Contains( ftree, ei->he->face, ei->he->flip->face ) &&
	  !Contains( vtree, ei->he->vertex, ei->he->flip->vertex ) ){
	// BINGO: neither tree contains this edge

	// trace chains to root
	std::vector<HalfEdgeCIter> vir; TraceToRoot( vir, ei->he->vertex, vtree );
	std::vector<HalfEdgeCIter> vjr; TraceToRoot( vjr, ei->he->flip->vertex, vtree );

	// remove common postfix (the paths may be merging before
	// reaching the root)
	while( !vir.empty() && !vjr.empty() && vir.back() == vjr.back() ){
	  vir.pop_back(); vjr.pop_back();
	}

	// now make final chain
	std::vector<EdgeCIter> g; g.push_back(ei);
	// and put them on the generator
	for( unsigned int i = 0; i < vir.size(); i++ ) g.push_back( vir[i]->edge );
	for( unsigned int j = vjr.size(); j > 0; j-- ) g.push_back( vjr[j-1]->edge );
	gens.push_back( g );
      }
    }
    return boundary;
  }

  void Mesh::ComputeRotation1Form( const unsigned int n, const bool constK ){

    const bool boundary = FindGenerators();
    cerr << gens.size() << " generators\n";
    const unsigned int en = edges.size(), vn = vertices.size(), fn = faces.size();
    const unsigned int twog = 2-fn+en-vn; assert( twog % 2 == 0 );

    std::map<EdgeCIter, unsigned int> ind;
    unsigned int i = 0;
    for( EdgeIter ei = edges.begin(); ei != edges.end(); ei++, i++ ){
      ind.insert( pair<const EdgeCIter, unsigned int>( ei, i ) );
    }

    SparseMatrix<Real> D( en, en );
    DenseMatrix<Real> rhs( en, 1 ), w( en, 1 );

    const double k = constK ? (((2.-double(twog))*2.*DDGConstants::PI)/double(fn)) : 0;
    i = 0; // row counter
    // sum conditions for every triangle; if there are no holes we
    // need to skip one triangle; if there are holes we use all
    // triangles
    bool skip1Triangle = !boundary;
    for( FaceCIter fi = faces.begin(); fi != faces.end(); fi++ ){
      if( fi->isBoundary() || skip1Triangle ){ skip1Triangle = false; continue; }
      const HalfEdgeIter eij = fi->he, ejk = fi->he->next, eki = fi->he->next->next;
      // dw = \tilde{\Omega}-\Omega
      D(i,ind[eij->edge]) = eij->sign();
      D(i,ind[ejk->edge]) = ejk->sign();
      D(i,ind[eki->edge]) = eki->sign();
      if( constK ) rhs(i) = k - fi->K();
      else rhs(i) = double(fi->sing)*2.*DDGConstants::PI/double(n) - fi->K();
      i++;
    }
    // closure condition for every vertex; skip first vertex
    for( VertexIter vi = ++vertices.begin(); vi != vertices.end(); vi++, i++ ){
      HalfEdgeIter hi = vi->he;
      do{
	D(i,ind[hi->edge]) = hi->sign() * hi->edge->cot();
      }while( ( hi = hi->flip->next ) != vi->he );
      rhs(i) = 0;
    }

    // add generators
    for( std::vector<std::vector<EdgeCIter> >::const_iterator gi = gens.begin(); gi != gens.end(); gi++, i++ ){
      for( std::vector<EdgeCIter>::const_iterator li = gi->begin(); li != gi->end(); li++ ){
	D(i,ind[*li]) = 1;
      }
      // holonomies around generators; these could be non-zero if
      // desired, thought of course they need to satisfy the index sum
      // condition; most natural choice most of the time is zero
      rhs(i) = 0;
    }
    assert( i == en );

    solveSquare( D, w, rhs );

    // load solution back
    i = 0; for( EdgeIter ei = edges.begin(); ei != edges.end(); ei++, i++ ) ei->w = w(i);
  }
} // namespace DDG
