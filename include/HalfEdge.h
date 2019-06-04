// -*- Mode: c++ -*-
// -----------------------------------------------------------------------------
// libDDG -- HalfEdge.h
// -----------------------------------------------------------------------------
//
// HalfEdge is used to define mesh connectivity.  (See the documentation for a
// more in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_HALFEDGE_H
#define DDG_HALFEDGE_H

#include "Vector.h"
#include "Types.h"
#include "Complex.h"
#include "Quaternion.h"

namespace DDG
{
  class HalfEdge
  {
  public:
    HalfEdgeIter next;
    // points to the next halfedge around the current face

    HalfEdgeIter flip;
    // points to the other halfedge associated with this edge

    VertexIter vertex;
    // points to the vertex at the "tail" of this halfedge

    EdgeIter edge;
    // points to the edge associated with this halfedge

    FaceIter face;
    // points to the face containing this halfedge

    bool onBoundary;
    // true if this halfedge is contained in a boundary
    // loop; false otherwise

    Vector texcoord;
    // texture coordinates associated with the triangle corner at the
    // "tail" of this halfedge

    // retrieves transport angle along this HalfEdge from parent
    // Edge depending on orientation
    double r( void ) const;

    // if we work with a changed connection \tilde{\nabla} = \nabla - Jw
    // we find the rotation form here
    double w( void ) const; // REMOVE w

    // sometimes we need to know the orientation of this halfedge with
    // respect to its parent edge
    int sign( void ) const;

    // retrieves mass matrix coefficient along this HalfEdge from
    // parent Edge depending on orientation
    Complex m( void ) const;

    // retrieves Energy matrix coefficient along this HalfEdge from
    // parent Edge depending on orientation
    Complex Es( void ) const;

    // accumulator function for mass matrix coefficient along this
    // HalfEdge updating parent Edge depending on orientation
    void mAcc( const Complex& m );

    // accumulator function for Energy matrix coefficient along this
    // HalfEdge updating parent Edge depending on orientation
    void EsAcc( const Complex& Es );

    // embedded geometry of this edge
    Vector geom( void ) const;

    // good old fashioned cotan of angle across from this halfedge
    const double cot( void ) const;

    // angle of triangle corner at the "tail" of this halfedge
    const double cornerAngle( void ) const;
  };
}

#endif

