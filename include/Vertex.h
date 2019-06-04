// -*- Mode: c++ -*-
// -----------------------------------------------------------------------------
// libDDG -- Vertex.h
// -----------------------------------------------------------------------------
//
// Vertex stores attributes associated with a mesh edge.  The iterator he
// points to its "outgoing" halfedge.  (See the documentation for a more
// in-depth discussion of the halfedge data structure.)
// 

#ifndef DDG_VERTEX_H
#define DDG_VERTEX_H

#include "Vector.h"
#include "Types.h"
#include "Complex.h"
#include "HalfEdge.h"

namespace DDG
{
  class Vertex
  {
  public:
    HalfEdgeIter he;
    // points to the "outgoing" halfedge
    HalfEdgeCIter parent;
    // needed in generators construction (shouldn't be here...)

    Vector position;
    // location of vertex in Euclidean 3-space

    // PS: K-vec/dir stuff
    // coefficient of PL section
    Complex u;

    // basis direction, which is a half edge
    const HalfEdgeIter X( void ) const { return he; }

    // resp. an actual vector in world space
    const Vector Xvector( void ) const;

    // angle scale factor for intrinsic tangent spaces
    double s;

    // mass of this vertex (happens to be real, but otherwise is a
    // complex coefficient); also called m_{ii} in the paper
    double m;

    // energy matrix coefficient Es_{ii}
    double Es;

    // Coefficient of PL section used in alignemnt (initialized to
    // Hopf differential as that is the canonical use of this)
    Complex q;

    // returns the kVec in world space
    const Vector kVec( void ) const;

    Vector normal;
    void updateNormal( void );
    // returns the vertex unit normal (volume gradient normal,
    // i.e., area weighted adding of incident face normals)

    bool isIsolated( void ) const;
    // returns true if the vertex is not contained in any face or edge; false otherwise

    // PS: counts number of edges, not faces (makes a difference at the boundary)
    const unsigned int valence( void ) const;
    // returns the number of incident faces

    double AngleOfEdge( const HalfEdgeCIter he ) const;
    // returns angle from X() direction to he in the flattened tangent space

    void EnforceBoundaryHalfEdgeConvention( void );
    // some code is greatly simplified if we may assume that the he
    // pointer of a boundary vertex points to the first boundary edge
    // in CCW order; enforce this here

    // PS: gross. but this is gotta be good enough for now (amortized cost is constant)
    bool onBoundary( void ) const {
      HalfEdgeIter hi = he;
      do{ if( hi->onBoundary ) return true; }while( ( hi = hi->flip->next ) != he );
      return false;
    }

    double BoundaryNormalAngle();
    // for a boundary vertex, return angle corresponding to the
    // (inward) vector orthogonal to the domain boundary

    Complex BoundaryValue( unsigned int n );
    // complex value corresponding to the boundary normal, for an n-direction field

    int id;
    // unique ID (used to index matrices)
  };
}

#endif

