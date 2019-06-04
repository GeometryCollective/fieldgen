#include <vector>
#include <math.h>
#include <cassert>
using namespace std;

#include "Vertex.h"
#include "Mesh.h"
#include "HalfEdge.h"
#include "Quaternion.h"

namespace DDG
{
  // returns the vertex normal
  void Vertex::updateNormal( void ) {
    Vector N(0,0,0);

    HalfEdgeCIter h = he;
      // PS: bugfix
    do{ if( !h->face->isBoundary() ) N += h->face->Anormal; }while( ( h = h->flip->next ) != he );

    // returns the volume gradient normal
    normal = N.unit();
  }

  vector<HalfEdge> isolated; // all isolated vertices point to isolated.begin()

  // returns true if the vertex is not contained in any face or edge; false otherwise
  bool Vertex::isIsolated( void ) const { return he == isolated.begin(); }

  // returns the number of incident faces
  const unsigned int Vertex :: valence( void ) const {
    unsigned int n = 0;

    HalfEdgeCIter h = he;
    do{ n++; }while( ( h = h->flip->next ) != he );
     
    return n;
  }

  // return the kVec as a worldspace vector
  const Vector Vertex::kVec( void ) const {
    const Vector x = Xvector();
    const Vector n = normal;
    // the corresponding tangent vector
    const Vector xt = ( x - dot(x,n)*n ).unit();
    const double theta = u.arg();
    const Quaternion r(cos(theta/2),sin(theta/2)*n);

    return u.norm()*(r*Quaternion(0,xt)*r.conj()).im();
  }

  const Vector Vertex::Xvector( void ) const{
    return X()->next->vertex->position - X()->vertex->position;
  }

  double
  Vertex::AngleOfEdge( const HalfEdgeCIter he ) const{
    // measure the angle from X() to he
    HalfEdgeIter hi = X(); double aCi = 0;
    while( hi != he ){
      // no angle to add for outer face; notice that we exit. The
      // assumption is that for a vertex on the boundary the
      // (only!) outer face will be the last face; so just exit
      if( hi->face->isBoundary() ) break;

      // get the angle that needs to be added to get to the next edge
      aCi += hi->face->alpha[hi->face->i(this->he->vertex)];

      // iterate
      hi = hi->next->next->flip;

      // guard against infinite loop
      assert( hi != X() );
    }
    return aCi * s;
  }

  void
  Vertex::EnforceBoundaryHalfEdgeConvention( void ){
    // if a vertex is on the boundary we want its he pointer to point
    // to the first boundary edge in CCW order; enforce this here
    HalfEdgeIter hi = he;
    do{
      // if there is a boundary face bend the he pointer and exit
      if( hi->flip->face->isBoundary() ){ he = hi; break; }
    }while( ( hi = hi->flip->next ) != he );
  }

  double Vertex::BoundaryNormalAngle()
  {
     HalfEdgeIter h = he;
     do
     {
        if( h->onBoundary ) break;
        h = h->next->next->flip;
     }
     while( h != he );

     return AngleOfEdge( h )/2.;
  }

  Complex Vertex::BoundaryValue( unsigned int n )
  {
     double theta = n * BoundaryNormalAngle();
     return Complex( cos(theta), sin(theta) );
  }
}

