// -----------------------------------------------------------------------------
// libDDG -- Viewer.h
// -----------------------------------------------------------------------------
//
// Viewer provides a graphical user interface (GUI) for inspecting and
// interacting with a Mesh object.  Viewer methods are static in order
// to make them compatible with GLUT callbacks.
//

#ifndef DDG_VIEWER_H
#define DDG_VIEWER_H

#include <string>
#include <GLUT/glut.h>
#include "Mesh.h"
#include "Camera.h"
#include "Shader.h"

namespace DDG
{
   class Viewer
   {
      public:
         static void init( void );
         // displays the viewer until the program ends

         static Mesh mesh;
         // surface mesh visualized by Viewer

      protected:
         // draw routines
         static void initGL( void );
         static void initGLUT( void );
         static void initGLSL( void );
         static void initLighting( void );
         static void drawSurface( void );
         static void drawMesh( void );
         static void setMeshMaterial( void );
         static void updateDisplayList( void );
	 static void drawField( void );
         static void drawPolygons( void );
         static void drawWireframe( void );
         static void drawHedgehog( void );
         static void drawIsolatedVertices( void );
         static void drawSingularities( void );
         static void drawInfo( void );

         enum stringAlignment
         {
            alignLeft,
            alignCenter,
            alignRight
         };
         static void drawString( std::string s, int x, int y, stringAlignment alignment = alignLeft );
         
         // GLUT callbacks
         static void display( void );
         static void idle( void );
         static void keyboard( unsigned char c, int x, int y );
         static void special( int i, int x, int y );
         static void mouse( int button, int state, int x, int y );
         static void motion( int x, int y );
         static void menu( int value );
         static void view( int value );
         
         // menu functions
	 static void mSmoothField( void );
	 static void mToggleAlignment( void );
	 static void mToggleFixedBoundary( void );
	 static void mToggleHedgehog( void );
	 static void mToggleNormalized( void );
         static void mResetMesh( void );
         static void mWriteMesh( void );
         static void mExit( void );
         static void mSmoothShaded( void );
         static void mTextured( void );
         static void mWireframe( void );
         static void mZoomIn( void );
         static void mZoomOut( void );
         static void mScreenshot( void );

         // unique identifiers for menus
         enum
         {
	    menuSmoothField,
	    menuInfConfDefo,
            menuResetMesh,
            menuWriteMesh,
            menuExit,
            menuSmoothShaded,
            menuWireframe,
            menuZoomIn,
            menuZoomOut,
            menuScreenshot,
            menuToggleAlignment,
            menuToggleFixedBoundary,
            menuToggleNormalized,
            menuToggleHedgehog
         };

         // draw state
         enum RenderMode
         {
            renderShaded = 0,
            renderWireframe = 1
         };

         static RenderMode mode;
         // current render mode
	 static bool fieldViz;
         static bool hedgehogViz;

         static void storeViewerState( void );
         static void restoreViewerState( void );
         static int windowSize[2];

         static Camera camera;
         // keeps track of view state

         static GLuint surfaceDL;
         // display list for mesh
         
         static Shader shader;
         // shader used to determine appearance of surface

         // field parameters
         static int fieldDegree;  // degree k of k-field
         static bool normalized;  // whether to normalize vectors
         static bool align;       // toggles alignment with curvature
         static bool fixBoundary; // toggles fixed boundary vectors (Dirichlet conditions)
         static double t;         // amount of alignment
         static double s;         // determines smoothness energy (between -1,1)
   };
}

#endif

