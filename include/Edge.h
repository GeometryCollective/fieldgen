// -*- Mode: c++ -*-
// -----------------------------------------------------------------------------
// libDDG -- Edge.h
// -----------------------------------------------------------------------------
//
// Edge stores attributes associated with a mesh edge.  The iterator he points
// to one of its two associated halfedges.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_EDGE_H
#define DDG_EDGE_H

#include "Types.h"
#include "Quaternion.h"
#include "Complex.h"

namespace DDG
{
  class Edge
  {
  public:
    // points to one of the two halfedges associated with this edge
    HalfEdgeIter he;

    // distinguished direction at this edge
    // PS: only needed for old code
    const HalfEdgeIter X( void ) const { return he; }
    const Vector Xvector( void ) const;

    // 1/2 sum of incident triangle cotans
    const double cot( void ) const;

    // transport along edge with orientation he
    double r;
    // additional \omega rotation form \tilde{\nabla} = \nabla - J\omega
    double w; // REMOVE w
    // mass matrix coefficient
    Complex m;
    // energy matrix coefficient
    Complex Es;
    // for Hopf alignment with orientation of he
    Complex q;
  };
}

#endif
