// -*- Mode: c++ -*-
// -----------------------------------------------------------------------------
// libDDG -- Face.h
// -----------------------------------------------------------------------------
//
// Face stores attributes associated with a mesh edge.  The iterator he points
// to one of its associated halfedges.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_FACE_H
#define DDG_FACE_H

#include "Types.h"
#include "Complex.h"
#include "Vector.h"

namespace DDG
{
  class Face
  {
  public:
    HalfEdgeIter he;
    // points to one of the halfedges associated with this face
    HalfEdgeCIter parent;
    // needed for generator construction (should not be here...)

    // PS: For a face the predicate is a function and called isBoundary
    // PS: for halfedges it is a member and called onBoundary
    // PS: I think making this consistent would be desirable
    bool isBoundary( void ) const;
    // returns true if this face corresponds to a
    // boundary loop; false otherwise

    // distinguished direction
    const HalfEdgeIter X( void ) const { return he; }
    const Vector Xvector( void ) const;

    // which is v? 0,1,2?
    unsigned int i( const VertexCIter v ) const;

    // field
    Complex u;

    // area of triangle
    double A;

    // corner angles
    double alpha[3];

    // index of this triangle
    int sing;

    Vector normal;
    void updateNormal( void );
    // returns the unit normal associated with this face; normal
    // orientation is determined by the circulation order of halfedges

    Vector Anormal;
    void updateANormal( void );
    // returns the area weighted normal

    double area( void ) const;
    // returns the triangle area

    Vector barycenter( void ) const;
    // returns the mean of vertex coordinates

    // as the name says
    void InitAreaAnglesAngleSums( void );

    // geometric curvature
    double K( void ) const;

    // PS: Wishlist
    // PS: unsigned int valence( void ) const; (or degree() if you wish)
    // PS: an assertion if things are not triangle in all the places that assume it

    // visualization ---------------------------------------------------
    Vector sampleUniform( void ) const;
    // samples a point uniformly at random, expressed in barycentric
    // coordinates relative to he->vertex, he->next->vertex, and
    // he->next->next->vertex, respectively.

    Vector sampleField( const Vector& p, double n ) const;
    // samples the interpolated field u at the point p expressed in
    // barycentric coordinates relative to he->vertex, he->next->vertex, and
    // he->next->next->vertex, respectively.
    // must specify the field degree n

    // converts from barycentric coordinates to (x,y,z) coords in R^3
    Vector toWorldCoordinates( const Vector& barycentric ) const;
  };
}

#endif
