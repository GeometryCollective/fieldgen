#include <algorithm>
#include <cmath>
#include <cassert>
#include <map>
#include <queue>
#include "Utility.h"
#include "Complex.h"
#include "Mesh.h"
#include "SparseMatrix.h"
#include "SectionIntegrals.h"

namespace DDG{

  // turn angle into unitary complex number
  inline const Complex Phase( const double a ) { return Complex(cos(a),sin(a)); }

  // angle from f to t in the plane orthogonal to N
  inline const double Angle( const Vector& f, const Vector& t, const Vector& N ){
    return atan2( dot(cross(f, t), N), dot(f, t) );
  }

  inline const double Clip( const double value, const double lower, const double upper ){
    return value < lower ? lower : ( value > upper ? upper : value );
  }

  void Mesh::ComputeConnectionAndHopf( void ){
    // at this point all face->alpha[] slots contain the Euclidean tip angles
    // now visit every edge and figure out transport along that edge
    for( EdgeIter ei = edges.begin(); ei != edges.end(); ei++ ){
      const VertexIter vi = ei->he->vertex, vj = ei->he->flip->vertex;
      
      // AngleOfEdge assumes that the alpha[] slots in all triangles
      // contain the Euclidean angles; it then returns the angle
      // measured from the distinguished direciton to the current edge
      // RESCALED (by s). So it is the flattened angle that is
      // measured here
      const double aCi = vi->AngleOfEdge( ei->he );
      // we subtract PI from this one since we don't want to know the
      // angle to ei->he->flip, but rather ei->he; undo the additional
      // PI rotation
      const double aCj = vj->AngleOfEdge( ei->he->flip ) - DDGConstants::PI;

      // what coefficient do I need to express X_i using the frame of X_j: X_i = r_{ij} X_j
      ei->r = ( aCj - aCi );
      // w holds possible modification of connection; initially: none
      ei->w = 0; // REMOVE w

      // let's compute q on this edge
      const double le = ei->he->geom().norm();
      if( ei->he->onBoundary || ei->he->flip->onBoundary ){
	// boundary edges have a zero dihedral angle
	ei->q = 0;
      }else{
	// need Clip() here to avoid being ever so slightly larger than 1 or smaller than -1
	// const double cs = Clip( dot(ei->he->face->normal, ei->he->flip->face->normal ), -1, 1 );
	// in the edge tangent space (distinguished direction is along edge) q is purely real
	// ei->q = -Complex(acos(cs)/2*le,0);
	ei->q = -Complex(Angle( ei->he->face->normal, ei->he->flip->face->normal, ei->he->geom().unit() )/2*le,0);
	// just making sure...
	assert( !isnan(ei->q.re) );
      }

      // assuming that the q is aligned with the edge it gets
      // transported into each end by the angles at the vertices with
      // the edge; since q is quadratic we use the square
      vj->q += Phase(aCj*2) * ei->q;
      vi->q += Phase(aCi*2) * ei->q;
    } // edge traversal
  }

  // new version which uses closed form integrals
  void Mesh::InitKVecDirData( void ){

    // we'll need transports along edges, areas (masses) of triangles
    // and masses (areas) of vertices; the latter is needed for the EV
    // problem as the mass matrix from the rhs gets incorporated into
    // the lhs so that we have just an ordinary EV problem
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      // initializations: mass, angle sums, Hopf differential
      vi->m = 0; vi->s = 0; vi->q = Complex(0,0);
    }

    for( FaceIter fi = faces.begin(); fi != faces.end(); fi++ ){
      if( fi->isBoundary() ) continue; // skip boundary face(s)
      // computes triangle area; puts Euclidean corner angles into
      // alpha[], adds them to the s accumulators in the three
      // vertices
      fi->InitAreaAnglesAngleSums();
    }
    
    // we used s to accumulate the total angle around a vertex; turn
    // it into the scale factor we need; on the boundary no flattening
    // is performed
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      vi->s = vi->onBoundary() ? 1 : (2*DDGConstants::PI / vi->s);
    }
    
    // assumes alpha angles and s scalings at vertices
    ComputeConnectionAndHopf();
  }

  void Mesh::ComputeEnergyAndMass( const unsigned int n, const double s ){

    // gotta zero out all working slots
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      vi->m = 0; vi->Es = 0;
    }
    for( EdgeIter ei = edges.begin(); ei != edges.end(); ei++ ){
      ei->m = Complex(0,0); ei->Es = Complex(0,0);
    }

    // each local stiffness matrix
    for( FaceIter fi = faces.begin(); fi != faces.end(); fi++ ){

      if( fi->isBoundary() ) continue; // nothing to do for boundary loop
      assert( fi->he->next->next->next == fi->he ); // better be a triangle

      // boundary of triangle ij, jk, ki
      const HalfEdgeIter he[3] = { fi->he, fi->he->next, fi->he->next->next };
      // get vertices for this triangle
      const VertexIter v[3] = { he[0]->vertex, he[1]->vertex, he[2]->vertex };
      // get edge differences, ij, jk, ki
      const Vector p[3] = { he[0]->geom(), he[1]->geom(), he[2]->geom() };
      // length squared
      const double pn2[3] = { p[0].norm2(), p[1].norm2(), p[2].norm2() };

      // transport coefficients along the edges (this is where
      // multiplication with k enters); we need to account for the
      // orientation of the halfedge relative to the edge; since we
      // only need the conjugate of these we'll conjugate right here
      const Complex rc[3] = // ADD w for flat bundle
	{ Phase(n*(he[0]->r()/*+he[0]->w()*/)).conj(),
	  Phase(n*(he[1]->r()/*+he[1]->w()*/)).conj(),
	  Phase(n*(he[2]->r()/*+he[2]->w()*/)).conj() };

      // area of this triangle
      const double A = fi->A;
      // line bundle curvature of this triangle
      const double om = n*(fi->K()/*+(he[0]->w()+he[1]->w()+he[2]->w())*/); // ADD w for flat line bundle

      // first the mass matrix
      const double Mii = A*MassII();
      // diagonal entries at vertices
      v[0]->m += Mii; v[1]->m += Mii; v[2]->m += Mii;

      const Complex Mij = A*MassIJ( om );
      // off diagonal at edges
      he[0]->mAcc( Mij*rc[0] ); he[1]->mAcc( Mij*rc[1] ); he[2]->mAcc( Mij*rc[2] );

      // energy matrix
      const double snKii = (s*(om/A))*Mii; assert( !isnan( snKii ) );
      // diagonal entries
      v[0]->Es += DirichletII( om, pn2[0], -dot(p[0],p[2]), pn2[2] )/A - snKii;
      v[1]->Es += DirichletII( om, pn2[1], -dot(p[1],p[0]), pn2[0] )/A - snKii;
      v[2]->Es += DirichletII( om, pn2[2], -dot(p[2],p[1]), pn2[1] )/A - snKii;

      const Complex snKij = s*((om/A)*Mij-Complex(0,.5));
      // off diagonal
      he[0]->EsAcc( (DirichletIJ( om, pn2[2], -dot(p[2],p[1]), pn2[1] )/A - snKij)*rc[0] );
      he[1]->EsAcc( (DirichletIJ( om, pn2[0], -dot(p[0],p[2]), pn2[2] )/A - snKij)*rc[1] );
      he[2]->EsAcc( (DirichletIJ( om, pn2[1], -dot(p[1],p[0]), pn2[0] )/A - snKij)*rc[2] );
    } // Face traversal
  }

  void Mesh::SetupEnergyMatrix( SparseMatrix<Complex> &A,
        			SparseMatrix<Complex> &M,
        			const unsigned int n, const double s,
                                double lambda )
  {
    const unsigned int nV = vertices.size();
    const unsigned int nE = edges.size();
    A.resize( nV, nV, nV+2*nE );
    M.resize( nV, nV, nV+2*nE );

    ComputeEnergyAndMass( n, s );
     
    // shift A by epsilon to avoid degeneracy
    // (does not affect smallest eigenvector)
    const double shift = -lambda + 1e-9;

    // temporary storage for matrix columns
    vector< vector<ColumnEntry> > cA( nV ), cM( nV );

    for( VertexIter vj  = vertices.begin(); vj != vertices.end(); ++vj )
    {
      size_t j = vj->id;

      vector<ColumnEntry>& Aj( cA[j] ); // nonzeros in jth column of A
      vector<ColumnEntry>& Mj( cM[j] ); // nonzeros in jth column of M
      Aj.reserve( 16 );
      Mj.reserve( 16 );

      Aj.push_back( ColumnEntry( j, Complex(vj->Es+shift*vj->m,0) ));
      Mj.push_back( ColumnEntry( j, Complex(vj->m,0) ));

      HalfEdgeCIter he = vj->he;
      do {
        EdgeCIter e = he->edge;
        VertexIter vi = he->flip->vertex;
        size_t i = vi->id;

        Complex Aij = e->Es + shift*e->m;
        Complex Mij = e->m;
        if( he->edge->he == he ){
          Aij = Aij.conj();
          Mij = Mij.conj();
        }

        Aj.push_back( ColumnEntry( i, Aij ));
        Mj.push_back( ColumnEntry( i, Mij ));

        he = he->flip->next;
      } while( he != vj->he );
    }

    for( unsigned int j = 0; j < nV; j++ )
    {
       A.appendCompressedColumn( cA[j] );
       M.appendCompressedColumn( cM[j] );
    }
  }

  void Mesh::SetupEnergyMatrixFixedBoundary( SparseMatrix<Complex> &A,
        			             SparseMatrix<Complex> &M,
        			             DenseMatrix<Complex> &b,
        			             const unsigned int n, const double s,
                                             double lambda )
  {
    const unsigned int nV = nInteriorVertices;
    const unsigned int nE = edges.size();
    A.resize( nV, nV, nV+2*nE );
    M.resize( nV, nV, nV+2*nE );
    b = DenseMatrix<Complex>( nV, 1 );
    b.zero( Complex(0.,0.) );

    ComputeEnergyAndMass( n, s );
     
    // shift A by epsilon to avoid degeneracy
    // (does not affect smallest eigenvector)
    const double shift = -lambda + 1e-9;

    // temporary storage for matrix columns
    vector< vector<ColumnEntry> > cA( nV ), cM( nV );

    for( VertexIter vj  = vertices.begin(); vj != vertices.end(); ++vj )
    {
       if( vj->onBoundary() ) continue;

       size_t j = vj->id;

       vector<ColumnEntry>& Aj( cA[j] ); // nonzeros in jth column of A
       vector<ColumnEntry>& Mj( cM[j] ); // nonzeros in jth column of M
       Aj.reserve( 16 );
       Mj.reserve( 16 );

       Aj.push_back( ColumnEntry( j, Complex(vj->Es+shift*vj->m,0) ));
       Mj.push_back( ColumnEntry( j, Complex(vj->m,0) ));

       HalfEdgeCIter he = vj->he;
       do {
          EdgeCIter e = he->edge;
          VertexIter vi = he->flip->vertex;

          Complex Aij = e->Es + shift*e->m;
          Complex Mij = e->m;
          if( he->edge->he == he ){
             Aij = Aij.conj();
             Mij = Mij.conj();
          }

          if( vi->onBoundary() )
          {
             // move boundary terms to the right-hand side
             b(j,0) -= Aij.conj() * vi->BoundaryValue(n);
          }
          else
          {
             size_t i = vi->id;
             Aj.push_back( ColumnEntry( i, Aij ));
             Mj.push_back( ColumnEntry( i, Mij ));
          }

          he = he->flip->next;
       } while( he != vj->he );
    }

    for( unsigned int j = 0; j < nV; j++ )
    {
       A.appendCompressedColumn( cA[j] );
       M.appendCompressedColumn( cM[j] );
    }
  }

  // energy minimizer
  void Mesh::ComputeSmoothest( const unsigned int n, const double s, const bool dir ){

    cerr << "Mesh::ComputeSmoothest: n: " << n << " s: " << s << " dir: " << dir << endl;

    double t0 = wallClock();

    const unsigned int nv = vertices.size();
    SparseMatrix<Complex> A( nv, nv ), M( nv, nv );

    SetupEnergyMatrix( A, M, n, s );

    // now find smallest eigenvector
    DenseMatrix<Complex> u(nv,1);
    u.randomize();

    smallestEigPositiveDefinite( A, M, u );
    cerr << "[eig] ev = " << rayleighQuotient( A, M, u ) << endl;

    // load result
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ )
    {
       vi->u = u(vi->id,0);
    }

    ComputeIndices( n );

    // take roots (and possibly normalize) vertex data
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      const Complex z = vi->u;
      vi->u = Phase(z.arg()/n) * ( dir ? 1. : pow(z.norm(),1./n) );
    }
    double t1 = wallClock(); printTiming( "compute smoothest field", t1-t0 );
    cerr << "triangles: " << faces.size() << endl;
  } // ComputeSmoothest

  void Mesh::ComputeSmoothestFixedBoundary( const unsigned int n, const double s, const bool dir ){

    cerr << "Mesh::ComputeSmoothestFixedBoundary: n: " << n << " s: " << s << " dir: " << dir << endl;

    double t0 = wallClock();

    const unsigned int nv = nInteriorVertices;
    SparseMatrix<Complex> A( nv, nv ), M( nv, nv );
    DenseMatrix<Complex> b;

    cout << "Build matrix" << endl;
    SetupEnergyMatrixFixedBoundary( A, M, b, n, s );

    cout << "Solve" << endl;
    // solve for smoothest field
    DenseMatrix<Complex> u(nv,1);
    solvePositiveDefinite( A, u, b );

    cerr << "Extract solution" << endl;
    // load result
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ )
    {
       if( !vi->onBoundary() ) // interior vertex
       {
          vi->u = u(vi->id,0);
       }
       else // boundary vertex
       {
          vi->u = vi->BoundaryValue(n);
       }
    }

    cout << "Compute indices" << endl;
    ComputeIndices( n );

    // take roots (and possibly normalize) vertex data
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      const Complex z = vi->u;
      vi->u = Phase(z.arg()/n) * ( dir ? 1. : pow(z.norm(),1./n) );
    }
    double t1 = wallClock(); printTiming( "compute smoothest field", t1-t0 );
    cerr << "triangles: " << faces.size() << endl;
  } // ComputeSmoothest

  // alignment with q slot at vertices code; typically q is the Hopf
  // differential, but could be something else; here we will assume
  // however that it is the Hopf differential
  double Mesh::SmoothestCurvatureAlignment( unsigned int n, double s, double lambda, bool dir ){
    cerr << "Mesh::CurvatureSmoothestAlignment: n: " << n << " s: " << s << " lambda: " << lambda << endl;
    if( n != 2 && n != 4 ){
      cerr << "alignment code requires n == 2 or n == 4; n is: " << n << endl;
      return 0;
    }
    double t0 = wallClock();

    const unsigned int nv = vertices.size();
    SparseMatrix<Complex> A( nv, nv ), M( nv, nv );

    SetupEnergyMatrix( A, M, n, s, lambda );

    // find solution to Poisson problem
    DenseMatrix<Complex> u(nv,1), q(nv,1);

    // load q; to simplify the t computation we need to normalize q
    unsigned int i = 0;
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); i++, vi++ ){
      q(vi->id,0) = ( n == 2 ? vi->q : vi->q*vi->q );
    }
    {
      u = M.multiply( q );
      double normQ = 0;
      for( i = 0; i < nv; i++ ) normQ += (q(i,0).conj()*u(i,0)).re;
      q = u/sqrt( normQ );
    }

    solvePositiveDefinite( A, u, q );

    // load result
    i = 0; for( VertexIter vi = vertices.begin(); vi != vertices.end(); i++, vi++ ) vi->u = u(vi->id,0);

    double normU = 0;
    {
      q = M*u;
      for( i = 0; i < nv; i++ ){
	normU += (u(i,0).conj()*q(i,0)).re;
      }
      normU = sqrt( normU );
    }

    ComputeIndices( n );

    // take roots (and possibly normalize) vertex data
    for( VertexIter vi = vertices.begin(); vi != vertices.end(); vi++ ){
      const Complex z = vi->u;
      vi->u = Phase(z.arg()/n) * ( dir ? 1. : pow(z.norm(),1./n) );
    }
    double t1 = wallClock(); printTiming( "compute aligned field", t1-t0 );
    cerr << "triangles: " << faces.size() << endl;

    // let the caller know what the actual t value was
    cerr << "corresponding t value: " << 1./(1.+normU) << endl;
    return 1./(1.+normU);
  }

  // determine singularity index of each triangle
  void Mesh::ComputeIndices( const unsigned int n ) {

    int totalIndex = 0;
    for( FaceIter fi = faces.begin(); fi != faces.end(); fi++ ){
      if( fi->isBoundary() ) continue;

      // boundary of triangle ij, jk, ki
      const HalfEdgeCIter he[3] = { fi->he, fi->he->next, fi->he->next->next };

      // index computation; first the geometric curvature of this
      // triangle
      const double sigma = fmodPI(n*fi->K()/*+n*(he[0]->w()+he[1]->w()+he[2]->w())*/); // ADD w for flat bundle

      // now walk around the boundary and look by how much the
      // vector field rotates.
      const Complex u[3] = { he[0]->vertex->u, he[1]->vertex->u, he[2]->vertex->u };

      // transport from vertex i to vertex j.
      const double rij[3] =  // ADD w for flat bundle
	{ fmodPI(n*(he[0]->r()/*+he[0]->w()*/)),
	  fmodPI(n*(he[1]->r()/*+he[1]->w()*/)),
	  fmodPI(n*(he[2]->r()/*+he[2]->w()*/)) };

      // now figure out how much the u field rotates all by itself;
      // transport a u_i to v_j and compare with u_j there
      const double aij[3] = {
	( u[0].unit() * Phase( rij[0] ) * u[1].unit().conj() ).arg(),
	( u[1].unit() * Phase( rij[1] ) * u[2].unit().conj() ).arg(),
	( u[2].unit() * Phase( rij[2] ) * u[0].unit().conj() ).arg() };

      //assert( abs( aij[0] + aij[1] + aij[2] - sigma ) < 4*DDGConstants::PI );

      // final result (note that we have an overall "sign error" since
      // we measured angles in CW orientation)
      const double s = -( aij[0] + aij[1] + aij[2] - sigma ) / (2*DDGConstants::PI );

      // this better be essentially an integer
      assert( abs( s - lround( s ) ) < 1e-6 );

      fi->sing = int( lround( s ) );

      totalIndex += fi->sing;

    } // Face traversal
  } // ComputeIndices

  void Mesh::clearSingularities( void )
  {
     for( FaceIter f = faces.begin(); f != faces.end(); f++ )
     {
        f->sing = 0;
     }
  }

} // namespace DDG
