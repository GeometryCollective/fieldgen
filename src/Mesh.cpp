#include <map>
#include <queue>
#include <fstream>
#include "Mesh.h"
#include "MeshIO.h"

using namespace std;

namespace DDG
{
   Mesh :: Mesh( void )
   {}
   
   Mesh :: Mesh( const Mesh& mesh )
   {
      *this = mesh;
   }
   
   class  HalfEdgeIterCompare { public: bool operator()( const  HalfEdgeIter& i, const  HalfEdgeIter& j ) const { return &*i < &*j; } };
   class HalfEdgeCIterCompare { public: bool operator()( const HalfEdgeCIter& i, const HalfEdgeCIter& j ) const { return &*i < &*j; } };
   class    VertexIterCompare { public: bool operator()( const    VertexIter& i, const    VertexIter& j ) const { return &*i < &*j; } };
   class   VertexCIterCompare { public: bool operator()( const   VertexCIter& i, const   VertexCIter& j ) const { return &*i < &*j; } };
   class      FaceIterCompare { public: bool operator()( const      FaceIter& i, const      FaceIter& j ) const { return &*i < &*j; } };
   class     FaceCIterCompare { public: bool operator()( const     FaceCIter& i, const     FaceCIter& j ) const { return &*i < &*j; } };
   class      EdgeIterCompare { public: bool operator()( const      EdgeIter& i, const      EdgeIter& j ) const { return &*i < &*j; } };
   class     EdgeCIterCompare { public: bool operator()( const     EdgeCIter& i, const     EdgeCIter& j ) const { return &*i < &*j; } };
   
   const Mesh& Mesh :: operator=( const Mesh& mesh )
   {
      map< HalfEdgeCIter, HalfEdgeIter, HalfEdgeCIterCompare > halfedgeOldToNew;
      map<   VertexCIter,   VertexIter,   VertexCIterCompare >   vertexOldToNew;
      map<     EdgeCIter,     EdgeIter,     EdgeCIterCompare >     edgeOldToNew;
      map<     FaceCIter,     FaceIter,     FaceCIterCompare >     faceOldToNew;
   
      // copy geometry from the original mesh and create a
      // map from pointers in the original mesh to
      // those in the new mesh
      halfedges.clear(); for( HalfEdgeCIter he = mesh.halfedges.begin(); he != mesh.halfedges.end(); he++ ) halfedgeOldToNew[ he ] = halfedges.insert( halfedges.end(), *he );
       vertices.clear(); for(   VertexCIter  v =  mesh.vertices.begin();  v !=  mesh.vertices.end();  v++ )   vertexOldToNew[ v  ] =  vertices.insert(  vertices.end(), *v  );
          edges.clear(); for(     EdgeCIter  e =     mesh.edges.begin();  e !=     mesh.edges.end();  e++ )     edgeOldToNew[ e  ] =     edges.insert(     edges.end(), *e  );
          faces.clear(); for(     FaceCIter  f =     mesh.faces.begin();  f !=     mesh.faces.end();  f++ )     faceOldToNew[ f  ] =     faces.insert(     faces.end(), *f  );
   
      // "search and replace" old pointers with new ones
      for( HalfEdgeIter he = halfedges.begin(); he != halfedges.end(); he++ )
      {
         he->next   = halfedgeOldToNew[ he->next   ];
         he->flip   = halfedgeOldToNew[ he->flip   ];
         he->vertex =   vertexOldToNew[ he->vertex ];
         he->edge   =     edgeOldToNew[ he->edge   ];
         he->face   =     faceOldToNew[ he->face   ];
      }
   
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ ) v->he = halfedgeOldToNew[ v->he ];
      for(   EdgeIter e =    edges.begin(); e !=    edges.end(); e++ ) e->he = halfedgeOldToNew[ e->he ];
      for(   FaceIter f =    faces.begin(); f !=    faces.end(); f++ ) f->he = halfedgeOldToNew[ f->he ];
   
      return *this;
   }

   int Mesh::read( const string& filename )
   {
      inputFilename = filename;
      ifstream in( filename.c_str() );

      if( !in.is_open() )
      {
         cerr << "Error reading from mesh file " << filename << endl;
         return 1;
      }

      int rval;
      if( !( rval = MeshIO::read( in, *this )))
      {
         indexVertices();
         updateGeometry();
      }
      return rval;
   }

   void Mesh::indexVertices()
   {
      int nVertices = 0;

      // index interior vertices first
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if( !v->onBoundary() )
         {
            v->id = nVertices;
            nVertices++;
         }
      }

      nInteriorVertices = nVertices;

      // then index boundary vertices
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         if( v->onBoundary() )
         {
            v->id = nVertices;
            nVertices++;
         }
      }
   }

   int Mesh::write( const string& filename ) const
   // reads a mesh from a Wavefront OBJ file; return value is nonzero
   // only if there was an error
   {
      ofstream out( filename.c_str() );

      if( !out.is_open() )
      {
         cerr << "Error writing to mesh file " << filename << endl;
         return 1;
      }

      MeshIO::write( out, *this );

      return 0;
   }

   bool Mesh::reload( void )
   {
      return read( inputFilename );
   }

   void Mesh::normalize( void )
   {
      // compute center of mass
      Vector c( 0., 0., 0. );
      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         c += v->position;
      }
      c /= (double) vertices.size();

      // translate to origin
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->position -= c;
      }

      // rescale such that the mesh sits inside the unit ball
      double rMax = 0.;
      for( VertexCIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         rMax = max( rMax, v->position.norm() );
      }
      for( VertexIter v = vertices.begin(); v != vertices.end(); v++ )
      {
         v->position /= rMax;
      }
   }

   // maximum distance from any vertex to the mesh centroid
   void Mesh :: updateRadius( void )
   {
      Vector c( 0., 0., 0. ); // centroid
      double A = 0.; // surface area

      // compute centroid as integral of position over surface area
      for( FaceCIter f  = faces.begin();
                     f != faces.end();
                     f ++ )
      {
	if( f->isBoundary() ) continue;
         double Af = f->Anormal.norm();
         Vector cf = f->barycenter();

         A += Af;
         c += Af*cf;
      }
      c /= A;

      // compute radius
      double r = 0.;
      for( VertexCIter v  = vertices.begin();
                       v != vertices.end();
                       v ++ )
      {
         r = max( r, (v->position-c).norm() );
      }

      radius = r;
   }

   void Mesh :: buildFaceTable( void )
   {
      // build alias table for probability distribution over faces
      for( size_t i = 0; i < faces.size(); i++ )
      {
         faceTable.push_back( faces[i].area() );
      }
      faceTable.build();
   }

   FaceCIter Mesh :: sampleUniform( void )
   {
      if( faceTable.size() == 0 )
      {
         buildFaceTable();
      }

      return faces[faceTable.sample()].he->face;
   }

   void Mesh :: updateNormals( void )
   {
      for( FaceIter f  = faces.begin();
                    f != faces.end();
                    f ++ )
      {
	if( f->isBoundary() ) continue;
         f->updateANormal();
         f->updateNormal();
      }

      for( VertexIter v  = vertices.begin();
                      v != vertices.end();
                      v ++ )
      {
         v->updateNormal();
      }
   }

   void Mesh :: updateGeometry( void )
   {
      // Note that order matters here, since
      // some routines use data cached in
      // previous routines.
      
      normalize();
      updateNormals();
      updateRadius();
   }
}

