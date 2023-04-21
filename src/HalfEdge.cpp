#include "HalfEdge.h"
#include "Edge.h"
#include "Mesh.h"
#include "Vertex.h"
#include "Vector.h"
#include <math.h>
#include <cassert>

namespace DDG
{
  // transport angle along this HalfEdge
  double HalfEdge::r( void ) const{
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    return edge->he == flip ? -edge->r : edge->r;
  }

  // additional rotation form // REMOVE w
  double HalfEdge::w( void ) const{
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    return edge->he == flip ? -edge->w : edge->w;
  }

  // orientation of this halfedge with respec to parent edge
  int HalfEdge::sign( void ) const{
    return edge->he == flip ? -1 : 1;
  }

  // mass matrix coefficient along this edge
  Complex HalfEdge::m( void ) const{
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    return edge->he == flip ? edge->m.conj() : edge->m;
  }

  // accumulator function for mass matrix coefficient along this edge
  void HalfEdge::mAcc( const Complex& m ){
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    if( edge->he == flip ) edge->m += m.conj(); else edge->m += m;
  }

  // energy matrix coefficient along this edge
  Complex HalfEdge::Es( void ) const{
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    return edge->he == flip ? edge->Es.conj() : edge->Es;
  }

  // accumulator function for Energy matrix coefficient along this edge
  void HalfEdge::EsAcc( const Complex& Es ){
    // if this halfedge is the one that the parent edge points to the
    // orientation is correct
    if( edge->he == flip ) edge->Es += Es.conj(); else edge->Es += Es;
  }

  Vector HalfEdge::geom( void ) const{
    return flip->vertex->position - vertex->position;
  }

  const double HalfEdge::cot( void ) const{
    // assumes triangles
    assert( face->he->next->next->next == face->he );

    const Vector pki = next->geom(), pij = next->next->geom();
    // <pij,pik>/|pij \times pik|
    return -dot(pij,pki)/cross(pki,pij).norm();
  }

  const double HalfEdge::cornerAngle( void ) const
  {
     Vector pi = vertex->position;
     Vector pj = next->vertex->position;
     Vector pk = next->next->vertex->position;
     Vector u = ( pj-pi ).unit();
     Vector v = ( pk-pi ).unit();
     return acos( dot( u, v ));
  }
}

