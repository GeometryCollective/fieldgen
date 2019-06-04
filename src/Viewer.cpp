#include <cmath>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
using namespace std;

#include "Viewer.h"
#include "Image.h"
#include "Utility.h"
#include "SectionIntegrals.h"

namespace DDG
{
   // declare static member variables
   Mesh Viewer::mesh;
   Viewer::RenderMode Viewer::mode = renderShaded;
   bool Viewer::fieldViz = false;
   bool Viewer::hedgehogViz = false;
   bool Viewer::showSingularities = true;
   GLuint Viewer::surfaceDL = 0;
   int Viewer::windowSize[2] = { 512, 512 };
   Camera Viewer::camera;
   Shader Viewer::shader;
   int Viewer::fieldDegree = 1;
   bool Viewer::normalized = true;
   bool Viewer::align = false;
   bool Viewer::fixBoundary = false;
   double Viewer::t = 0.;
   double Viewer::s = 0.;
   
   void Viewer :: init( void )
   {
      restoreViewerState();
      initGLUT();
      initGL();
      initGLSL();
   
      updateDisplayList();
   
      glutMainLoop();
   }

   void Viewer :: initGL( void )
   {
      glClearColor( 1., 1., 1., 0. );
      initLighting();
   }
   
   void Viewer :: initGLUT( void )
   {
      int argc = 0;
      vector< vector<char> > argv(1);
   
      // initialize window
      glutInitWindowSize( windowSize[0], windowSize[1] );
      glutInitDisplayString("rgba stencil double samples=16");
      //glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH );
      glutInit( &argc, (char**)&argv );
      glutCreateWindow( "DDG" );
   
      // specify callbacks
      glutDisplayFunc  ( Viewer::display  );
      glutIdleFunc     ( Viewer::idle     );
      glutKeyboardFunc ( Viewer::keyboard );
      glutSpecialFunc  ( Viewer::special  );
      glutMouseFunc    ( Viewer::mouse    );
      glutMotionFunc   ( Viewer::motion   );
   
      // initialize menus
      int viewMenu = glutCreateMenu( Viewer::view );
      glutSetMenu( viewMenu );
      glutAddMenuEntry( "[m] Smooth Shaded",         menuSmoothShaded        );
      glutAddMenuEntry( "[f] Wireframe",             menuWireframe           );
      glutAddMenuEntry( "[↑] Zoom In",               menuZoomIn              );
      glutAddMenuEntry( "[↓] Zoom Out",              menuZoomOut             );
      glutAddMenuEntry( "[*] Toggle Singularities",  menuToggleSingularities );

      int mainMenu = glutCreateMenu( Viewer::menu );
      glutSetMenu( mainMenu );
      glutAddMenuEntry( "[u] Smooth Field",          menuSmoothField         );
      glutAddMenuEntry( "[c] Toggle Curvature",      menuToggleAlignment     );
      glutAddMenuEntry( "[b] Toggle Fixed Boundary", menuToggleFixedBoundary );
      glutAddMenuEntry( "[n] Toggle Normalization",  menuToggleNormalized    );
      glutAddMenuEntry( "[r] Reset Mesh",            menuResetMesh           );
      glutAddMenuEntry( "[w] Write Mesh",            menuWriteMesh           );
      glutAddMenuEntry( "[\\] Screenshot",           menuScreenshot          );
      glutAddMenuEntry( "[esc] Exit",                menuExit                );
      glutAddSubMenu( "View", viewMenu );
      glutAttachMenu( GLUT_RIGHT_BUTTON );
   }

   void Viewer :: initGLSL( void )
   {
      shader.loadVertex( "shaders/vertex.glsl" );
      shader.loadFragment( "shaders/fragment.glsl" );
   }
   
   void Viewer :: initLighting( void )
   {
      GLfloat position[4] = { 20., 30., 40., 0. };
      glLightfv( GL_LIGHT0, GL_POSITION, position );
      glEnable( GL_LIGHT0 );
      glEnable( GL_NORMALIZE );
   }
   
   void Viewer :: menu( int value )
   {
      switch( value )
      {
         case( menuSmoothField ):
	    mSmoothField();
	    break;
         case( menuToggleHedgehog ):
            mToggleHedgehog();
            break;
         case( menuResetMesh ):
            mResetMesh();
            break;
         case( menuToggleAlignment ):
            mToggleAlignment();
            break;
         case( menuToggleFixedBoundary ):
            mToggleFixedBoundary();
            break;
         case( menuToggleNormalized ):
            mToggleNormalized();
            break;
         case( menuWriteMesh ):
            mWriteMesh();
            break;
         case( menuScreenshot ):
            mScreenshot();
            break;
         case( menuExit ):
            mExit();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: view( int value )
   {
      switch( value )
      {
         case( menuSmoothShaded ):
            mSmoothShaded();
            break;
         case( menuWireframe ):
            mWireframe();
            break;
         case( menuZoomIn ):
            mZoomIn();
            break;
         case( menuZoomOut ):
            mZoomOut();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: mSmoothField( void )
   {
      static bool first_call = true;
      if( first_call )
      {
         mesh.InitKVecDirData();
         first_call = false;
      }

      mesh.clearSingularities();

      if( align )
      {
         mesh.SmoothestCurvatureAlignment( fieldDegree, s, t, normalized );
      }
      else if( fixBoundary )
      {
         mesh.ComputeSmoothestFixedBoundary( fieldDegree, s, normalized );
      }
      else
      {
         mesh.ComputeSmoothest( fieldDegree, s, normalized );
      }

      fieldViz = true;

      updateDisplayList();
   }

   void Viewer :: mResetMesh( void )
   {
      mesh.reload();
      updateDisplayList();
   }
   
   void Viewer :: mWriteMesh( void )
   {
      mesh.write( "out.obj" );
   }
   
   void Viewer :: mExit( void )
   {
      storeViewerState();
      exit( 0 );
   }

   void Viewer :: mToggleAlignment( void )
   {
      align = !align;
      if( align ) fixBoundary = false;
   }
   
   void Viewer :: mToggleFixedBoundary( void )
   {
      fixBoundary = !fixBoundary;
      if( fixBoundary ) align = false;
   }
   
   void Viewer :: mToggleNormalized( void )
   {
      normalized = !normalized;
   }
   
   void Viewer :: mSmoothShaded( void )
   {
      mode = renderShaded;
      updateDisplayList();
   }
   
   void Viewer :: mWireframe( void )
   {
      mode = renderWireframe;
      updateDisplayList();
   }

   void Viewer :: mToggleSingularities( void )
   {
      showSingularities = !showSingularities;
      updateDisplayList();
   }

   void Viewer :: mZoomIn( void )
   {
      camera.zoomIn();
   }

   void Viewer :: mZoomOut( void )
   {
      camera.zoomOut();
   }

   void Viewer :: mScreenshot( void )
   {
      static int index = 0;
   
      // get window width and height
      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      int w = view[2];
      int h = view[3];
   
      // get pixels
      Image image( w, h );
      glReadPixels( 0, 0, w, h, GL_BGR, GL_FLOAT, &image(0,0) );
   
      stringstream filename;
      filename << "frames/viewer" << setw(8) << setfill( '0' ) << index << ".tga";
      image.write( filename.str().c_str() );
   
      index++;
   }
   
   void Viewer :: keyboard( unsigned char c, int x, int y )
   {
      switch( c )
      {
         case 'K':
            fieldDegree--;
            if( fieldDegree < 1 ) fieldDegree = 1;
            break;
         case 'k':
            fieldDegree++;
            break;
         case 's':
            s += .01;
            s = max( -1.+1e-7, min( 1.-1e-7, s ));
            break;
         case 'S':
            s -= .01;
            s = max( -1.+1e-7, min( 1.-1e-7, s ));
            break;
         case 't':
            t += .01;
            break;
         case 'T':
            t -= .01;
            break;
         case 'm':
            mSmoothShaded();
            break;
         case 'f':
            mWireframe();
            break;
         case 'w':
            mWriteMesh();
            break;
         case 'r':
            mResetMesh();
            break;
         case '\\':
            mScreenshot();
            break;
         case ' ':
	    mSmoothField();
	    break;
         case 'b':
	    mToggleFixedBoundary();
	    break;
         case 'c':
	    mToggleAlignment();
	    break;
         case '*':
	    mToggleSingularities();
	    break;
         // case 'n':
         //    mToggleNormalized();
         //    break;
         case 27:
            mExit();
            break;
         default:
            break;
      }

      // glutPostRedisplay();
   }

   void Viewer :: special( int i, int x, int y )
   {
      switch( i )
      {
         case GLUT_KEY_UP:
            camera.zoomIn();
            break;
         case GLUT_KEY_DOWN:
            camera.zoomOut();
            break;
         case 27:
            mExit();
            break;
         default:
            break;
      }
   }
   
   void Viewer :: display( void )
   {
      glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

      shader.enable();
   
      glMatrixMode( GL_PROJECTION );
      glLoadIdentity();
   
      GLint viewport[4];
      glGetIntegerv( GL_VIEWPORT, viewport );
      double aspect = (double) viewport[2] / (double) viewport[3];
      const double fovy = 50.;
      const double clipNear = .01;
      const double clipFar = 1000.;
      gluPerspective( fovy, aspect, clipNear, clipFar );
   
      glMatrixMode( GL_MODELVIEW );
      glLoadIdentity();

      Quaternion    eye = Vector( 0., 0., -2.5*camera.zoom );
      Quaternion center = Vector( 0., 0., 0. );
      Quaternion     up = Vector( 0., 1., 0. );

      gluLookAt(    eye[1],    eye[2],    eye[3],
                 center[1], center[2], center[3],
                     up[1],     up[2],     up[3] );


      Quaternion r = camera.currentRotation();

      eye = r.conj() * eye * r;
      GLint uniformEye = glGetUniformLocation( shader, "eye" );
      glUniform3f( uniformEye, eye[1], eye[2], eye[3] );

      Quaternion light = Vector( -1., 1., -2. );
      light = r.conj() * light * r;
      GLint uniformLight = glGetUniformLocation( shader, "light" );
      glUniform3f( uniformLight, light[1], light[2], light[3] );

      camera.setView();
   
      drawSurface();

      shader.disable();

      drawInfo();
   
      glutSwapBuffers();
   }

   void Viewer :: setMeshMaterial( void )
   {
      GLfloat  diffuse[4] = { 1., .5, .25, 1. };
      GLfloat specular[4] = { .3, .3, .30, 1. };
      GLfloat  ambient[4] = { .2, .2, .50, 1. };
   
      glColor4fv( diffuse );
      glMaterialfv( GL_FRONT_AND_BACK, GL_DIFFUSE,   diffuse  );
      glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR,  specular );
      glMaterialfv( GL_FRONT_AND_BACK, GL_AMBIENT,   ambient  );
      glMaterialf ( GL_FRONT_AND_BACK, GL_SHININESS, 16.      );
   }
   
   void Viewer :: drawSurface( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_DEPTH_TEST );
      glDisable( GL_LIGHTING );
   
      glCallList( surfaceDL );
   
      glPopAttrib();
   }
   
   void Viewer :: drawMesh( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glEnable( GL_POLYGON_OFFSET_FILL );
      glPolygonOffset( 1., 1. );
   
      drawPolygons();
   
      glDisable( GL_POLYGON_OFFSET_FILL );
   
      if( mode == renderWireframe )
      {
         shader.disable();
         drawWireframe();
	 // shader.enable();
      }

      //if( hedgehogViz )
      //{
      //   shader.disable();
      //   drawHedgehog();
      //}

      drawIsolatedVertices();

      glPopAttrib();
   }

   void Viewer :: drawPolygons( void )
   {
      glBegin( GL_TRIANGLES );
      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         if( f->isBoundary() ) continue;

         if( mode == renderWireframe )
         {
            Vector N = f->normal;
            glNormal3dv( &N[0] );
         }

         HalfEdgeCIter he = f->he;
         do
         {
            if( mode != renderWireframe )
            {
               Vector N = he->vertex->normal;
               glNormal3dv( &N[0] );
            }

            glVertex3dv( &he->vertex->position[0] );

            he = he->next;
         }
         while( he != f->he );
      }
      glEnd();
   }

   void Viewer :: drawField( void )
   {
      const double scale = (fieldDegree==1) ? 1. : .5;

      glPushAttrib( GL_ALL_ATTRIB_BITS );

      shader.disable();

      glLineWidth(1.);
      glColor4f( 0., 0., 0., .5 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

      for( VertexIter vi = mesh.vertices.begin(); vi != mesh.vertices.end(); vi++ ){
	const Vector x = vi->Xvector(), c = vi->position, n = vi->normal;
	const Vector rx = x.norm() * cross(n,cross(x,n)).unit();
	const double theta = vi->u.arg();
	for( int i = 0; i < fieldDegree; i++ ){
	   const Quaternion r(cos(theta/2+i*DDGConstants::PI/fieldDegree),sin(theta/2+i*DDGConstants::PI/fieldDegree)*n);
	   const Vector L = scale*vi->u.norm()*(r*Quaternion(0,rx)*r.conj()).im();
	   glBegin( GL_LINES );
	   glVertex3dv( &c[0] );
	   glVertex3dv( &(c+L)[0] );
	   glEnd();
	 }
      }

      glPopAttrib();
   }

   void Viewer :: drawWireframe( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      glColor4f( 0., 0., 0., .5 );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );

      glLineWidth( 10.0f );
      glBegin( GL_LINES );
      for( std::vector< std::vector<EdgeCIter> >::const_iterator el = mesh.gens.begin(); el != mesh.gens.end(); el++ ){
	for( std::vector<EdgeCIter>::const_iterator e = el->begin(); e != el->end(); e++ ){
	  glVertex3dv( &((*e)->he->vertex->position[0]) );
	  glVertex3dv( &((*e)->he->flip->vertex->position[0]) );
	}
      }
      glEnd();
      glLineWidth( 1.0f );
      glBegin( GL_LINES );
      for( EdgeCIter e  = mesh.edges.begin();
            e != mesh.edges.end();
            e ++ )
      {
         glVertex3dv( &e->he->vertex->position[0] );
         glVertex3dv( &e->he->flip->vertex->position[0] );
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: drawIsolatedVertices( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );

      // draw with big, round, red dots
      glPointSize( 5 );
      glHint( GL_POINT_SMOOTH_HINT, GL_NICEST );
      glEnable( GL_POINT_SMOOTH );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glColor4f( 1., 0., 0., 1. ); // red

      glBegin( GL_POINTS );
      for( VertexCIter v  = mesh.vertices.begin();
                       v != mesh.vertices.end();
                       v ++ )
      {
         if( v->isIsolated() )
         {
            glVertex3dv( &v->position[0] );
         }
      }
      glEnd();

      glPopAttrib();
   }

   void Viewer :: drawSingularities( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );
      glMatrixMode( GL_MODELVIEW );

      shader.enable();

      double radius = .005 * mesh.radius;
      int nSlices = 24;
      int nStacks = 24;

      for( FaceCIter f  = mesh.faces.begin();
                     f != mesh.faces.end();
                     f ++ )
      {
         if( f->isBoundary() ) continue;
         if( f->sing != 0 )
         {
            Vector c = f->barycenter();

            if( f->sing < 0 ) glColor3f( 0., 0., 1. );
            else              glColor3f( 1., 0., 0. );

            glPushMatrix();
            glTranslatef( c.x, c.y, c.z );
            glutSolidSphere( radius, nSlices, nStacks );
            glPopMatrix();
         }
      }

      glPopAttrib();
   }
   
   void Viewer :: updateDisplayList( void )
   {
      if( surfaceDL )
      {
         glDeleteLists( surfaceDL, 1 );
         surfaceDL = 0;
      }
   
      surfaceDL = glGenLists( 1 );
   
      setMeshMaterial();
   
      glNewList( surfaceDL, GL_COMPILE );
      drawMesh();
      if( fieldViz && !hedgehogViz )
      {
         drawField();
      }
      if( fieldViz || hedgehogViz )
      {
         if( showSingularities )
         {
            drawSingularities();
         }
      }
      glEndList();

      display();
   }
   
   void Viewer :: mouse( int button, int state, int x, int y )
   {
      camera.mouse( button, state, x, y );
   }

   void Viewer :: motion( int x, int y )
   {
      camera.motion( x, y );
   }
   
   void Viewer :: idle( void )
   {
      camera.idle();
      glutPostRedisplay();
   }

   void Viewer :: storeViewerState( void )
   {
      ofstream out( ".viewer_state.txt" );

      out << camera.rLast[0] << endl;
      out << camera.rLast[1] << endl;
      out << camera.rLast[2] << endl;
      out << camera.rLast[3] << endl;

      GLint view[4];
      glGetIntegerv( GL_VIEWPORT, view );
      out << view[2] << endl;
      out << view[3] << endl;

      out << (int) mode << endl;
   }

   void Viewer :: restoreViewerState( void )
   {
      ifstream in( ".viewer_state.txt" );

      if( !in.is_open() ) return;

      in >> camera.rLast[0];
      in >> camera.rLast[1];
      in >> camera.rLast[2];
      in >> camera.rLast[3];

      in >> windowSize[0];
      in >> windowSize[1];

      int m;
      in >> m;
      mode = (RenderMode) m;
   }

   void Viewer :: drawString( string s, int x, int y, stringAlignment alignment )
   {
      if( alignment == alignRight )
      {
         glRasterPos2f( x-9*(s.length()+1), y );
      }
      else if( alignment == alignCenter )
      {
         glRasterPos2f( x-9*(s.length()/2+1), y );
      }
      else
      {
         glRasterPos2f( x, y );
      }

      for( size_t i = 0; i < s.length(); i++ )
      {
         glutBitmapCharacter( GLUT_BITMAP_9_BY_15, s[i] );
      }
   }

   void Viewer :: drawInfo( void )
   {
      glPushAttrib( GL_ALL_ATTRIB_BITS );
      glDisable( GL_DEPTH_TEST );

      GLint viewport[4];
      glGetIntegerv( GL_VIEWPORT, viewport );
      int W = viewport[2];
      int H = viewport[3];

      glMatrixMode( GL_PROJECTION ); glPushMatrix(); glLoadIdentity();
      glOrtho( 0., W, 0., H, .1, 100. );
      glMatrixMode( GL_MODELVIEW  ); glPushMatrix(); glLoadIdentity();
      glTranslatef( 0., 0., -1. );

      glColor3f( 0., 0., 0. );

      int h = 16;
      int hInc = 14;
      
      // display field degree
      {
         stringstream ss;
         ss << "k: " << fieldDegree;
         drawString( ss.str(), 16, H-h, alignLeft );
         h += hInc;
      }
      
      // display curvature alignment
      {
         stringstream ss;
         ss << "curvature alignment: " << (align?"on":"off");
         drawString( ss.str(), 16, H-h );
         h += hInc;
      }
      
      // display boundary alignment
      {
         stringstream ss;
         ss << "boundary alignment: " << (fixBoundary?"on":"off");
         drawString( ss.str(), 16, H-h );
         h += hInc;
      }
      
      // display field degree
      {
         stringstream ss;
         ss << "normalized: " << (normalized?"true":"false");
         drawString( ss.str(), 16, H-h );
         h += hInc;
      }

      // display s
      {
         stringstream ss;
         ss << "s: " << s;
         drawString( ss.str(), 16, H-h );
         h += hInc;
      }

      // display t
      {
         stringstream ss;
         ss << "t: " << t;
         drawString( ss.str(), 16, H-h );
         h += hInc;
      }

      glMatrixMode( GL_PROJECTION ); glPopMatrix();
      glMatrixMode( GL_MODELVIEW  ); glPopMatrix();

      glPopAttrib();
   }

   void Viewer :: mToggleHedgehog( void )
   {
      hedgehogViz = !hedgehogViz;
      updateDisplayList();
   }

   void Viewer :: drawHedgehog( void )
   {
      const int nSamples = 10000;
      const double arrowLength = .05;

      glPushAttrib( GL_ALL_ATTRIB_BITS );
      glEnable( GL_BLEND );
      glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
      glLineWidth( 2. );

      glBegin( GL_LINES );
      for( int sample = 0; sample < nSamples; sample++ )
      {
         FaceCIter f = mesh.sampleUniform();
         Vector p = f->sampleUniform();
         Vector Z = f->sampleField( p, fieldDegree );
         Vector a = f->toWorldCoordinates( p );
         Vector b = a + Z*arrowLength;

         glColor4f( 0., 0., 0., 1. ); glVertex3dv( &a[0] );
         glColor4f( 0., 0., 0., 0. ); glVertex3dv( &b[0] );
      }
      glEnd();
      
      glPopAttrib();
   }
}

