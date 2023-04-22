// -----------------------------------------------------------------------------
// libDDG -- MeshIO.h
// -----------------------------------------------------------------------------
//
// MeshIO handles input/output operations for Mesh objects.  Currently the only
// supported mesh format is Wavefront OBJ -- for a format specification see
//
//   http://en.wikipedia.org/wiki/Wavefront_.obj_file
//
// Note that vertex normals and material properties are currently ignored.
//

#ifndef DDG_MESHIO_H
#define DDG_MESHIO_H

#include <iosfwd>
#include <string>
#include <sstream>
#include <vector>

#include "Vector.h"

namespace DDG
{
   class Mesh;
   class Index;
   class MeshData;

   class MeshIO
   {
      public:
         static int read( std::istream& in, Mesh& mesh );
         // reads a mesh from a valid, open input stream in

         static void write( std::ostream& out, const Mesh& mesh, unsigned int n );
         // writes a mesh to a valid, open output stream out

         static  int buildMesh( const MeshData& data, Mesh& mesh );
	 // made visible to help interface with Houdini
      protected:
         static  int readMeshData( std::istream& in, MeshData& data );
         static void readPosition( std::stringstream& ss, MeshData& data );
         static void readTexCoord( std::stringstream& ss, MeshData& data );
         static void readNormal  ( std::stringstream& ss, MeshData& data );
         static void readFace    ( std::stringstream& ss, MeshData& data );
         static void readAlignment( std::stringstream& ss, MeshData& data );
         static void readFaceAlignment( std::stringstream& ss, MeshData& data );
         static Index parseFaceIndex( const std::string& token );
         static void preallocateMeshElements( const MeshData& data, Mesh& mesh );
         static void checkIsolatedVertices( const Mesh& Mesh );
         static void checkNonManifoldVertices( const Mesh& Mesh );
   };

   class Index
   {
      public:
         Index( void )
         {}
   
         Index( int p, int t, int n )
         : position( p ), texcoord( t ), normal( n )
         {}
   
         bool operator<( const Index& i ) const
         {
            if( position < i.position ) return true;
            if( position > i.position ) return false;
            if( texcoord < i.texcoord ) return true;
            if( texcoord > i.texcoord ) return false;
            if(   normal < i.normal   ) return true;
            if(   normal > i.normal   ) return false;
            return false;
         }
   
         int position;
         int texcoord;
         int normal;
   };
   
   class MeshData
   {
      public:
         std::vector<Vector> positions;
         std::vector<Vector> texcoords;
         std::vector<Vector> normals;
         std::vector<Vector> alignments;
         std::vector< std::vector< Index > > indices;
   };
}

#endif

