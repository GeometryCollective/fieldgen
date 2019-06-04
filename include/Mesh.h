// -----------------------------------------------------------------------------
// libDDG -- Mesh.h
// -----------------------------------------------------------------------------
//
// Mesh represents a polygonal surface mesh using the halfedge data structure.
// It is essentially a large collection of disjoint vertices, edges, and faces
// that are ``glued together'' by halfedges which encode connectivity (see
// the documentation for an illustration).  By construction, the halfedge data
// structure cannot represent nonorientable surfaces or meshes with nonmanifold
// edges.
//
// Mesh elements are referenced using iterators -- common usage of these
// iterators is to either traverse an entire vector of mesh elements:
//
//    // visit all vertices
//    for( VertexIter i = vertices.begin(); i != vertices.end(); i++ )
//    {
//       //...
//    }
//
// or to perform a local traversal over the neighborhood of some mesh element:
//
//    // visit both halfedges of edge e
//    HalfEdgeIter he = e->he;
//    do
//    {
//       // ...
//
//       he = he->flip;
//    }
//    while( he != e->he );
//
// (See Types.h for an explicit definition of iterator types.)
//
// Meshes with boundary are handled by creating an additional face for each
// boundary loop (the method Face::isBoundary() determines whether a given
// face is a boundary loop).  Isolated vertices (i.e., vertiecs not contained
// in any edge or face) reference a dummy halfedge and can be checked via
// the method Vertex::isIsolated().
//

#ifndef DDG_MESH_H
#define DDG_MESH_H

#include <vector>
#include <string>

#include "HalfEdge.h"
#include "Vertex.h"
#include "Edge.h"
#include "Face.h"
#include "AliasTable.h"

namespace DDG
{
  typedef std::vector<HalfEdgeCIter> Chain;

  class Mesh
  {
  public:
    Mesh( void );
    // constructs an empty mesh

    Mesh( const Mesh& mesh );
    // constructs a copy of mesh

    const Mesh& operator=( const Mesh& mesh );
    // copies mesh

    int read( const std::string& filename );
    // reads a mesh from a Wavefront OBJ file; return value is nonzero
    // only if there was an error

    int write( const std::string& filename, unsigned int n ) const;
    // writes a mesh to a Wavefront OBJ file; return value is nonzero
    // only if there was an error
    // the value n gives the degree of the field

    bool reload( void );
    // reloads a mesh from disk using the most recent input filename
         
    void normalize( void );
    // centers around the origin and rescales to have unit radius

    void updateGeometry( void );
    // run updateGeometry() after changing vertex positions

    void findGenerators( void );
    void buildDualSpanningTree( void );
    void buildPrimalCoTree( void );

    // PS: k-vec/dir stuff
    // initializes up all manner of needed slots
    void InitKVecDirData( void ); 
    // solves for the smoothest of degree n, energy parameter s,
    // direction? (or vector)
    void ComputeSmoothest( const unsigned int n, const double s, const bool dir );
    // solves for smoothest field with boundary vectors fixed to boundary normals
    void ComputeSmoothestFixedBoundary( const unsigned int n, const double s, const bool dir );
    // solver function when using alignment energy
    double SmoothestCurvatureAlignment( const unsigned int n, const double s,
					const double lambda, const bool dir );

    // computes energy and mass matrices into the mesh to then be
    // loaded later into the global matrix; gets called by
    // SetupEnergyMatrix;
    void ComputeEnergyAndMass( const unsigned int n, const double s );
    // loads the matrix (gets called by ComputeSmoothest)
    void SetupEnergyMatrix( SparseMatrix<DDG::Complex> &A, SparseMatrix<DDG::Complex> &M,
			    const unsigned int n, const double s, double lambda = 0. );
    // loads the matrix (gets called by ComputeSmoothestFixedBoundary)
    void SetupEnergyMatrixFixedBoundary( SparseMatrix<DDG::Complex> &A, SparseMatrix<DDG::Complex> &M, DenseMatrix<DDG::Complex> &b,
			    const unsigned int n, const double s, double lambda = 0. );

    // computes indices of all triangles (gets called by ComputeSmoothest)
    void ComputeIndices( const unsigned int n );

    // uses q slot to do alignment (gets called by InitKVecDirData)
    void ComputeConnectionAndHopf( void );

    std::vector<HalfEdge> halfedges;
    std::vector<Vertex>   vertices;
    std::vector<Edge>     edges;
    std::vector<Face>     faces;
    std::vector< std::vector<EdgeCIter> > gens;
    // storage for mesh elements

    // returns maximum distance from any vertex to the mesh centroid
    double radius;

    FaceCIter sampleUniform( void );
    // returns a face uniformly at random relative to the surface area measure

    void clearSingularities( void );

  protected:
    std::string inputFilename;

    int nInteriorVertices;
    void indexVertices();
    // assign unique id to each vertex (interior first, then boundary)
    
    void updateNormals( void );
    void updateRadius( void );

    void buildFaceTable( void );
    AliasTable faceTable;
    // used to sample faces uniformly at random (for visualization)
  };
}

#endif

