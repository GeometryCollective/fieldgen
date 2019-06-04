#include "Face.h"
#include "Mesh.h"
#include "Vector.h"
#include "Utility.h"
#include <cassert>
#include <cmath>
#include <iostream>

using namespace std;

namespace DDG
{
  // unit normal
  void Face::updateNormal( void )
  {
     normal = Anormal.unit();
  }

  // area weighted normal
  void Face::updateANormal( void )
  {
    assert( he->next->next->next == he ); // assuming triangles

    const Vector p0 = he->vertex->position,
      p1 = he->next->vertex->position,
      p2 = he->next->next->vertex->position;

    Anormal = cross( p1-p0, p2-p0 )/2;
  }

  double Face::area( void ) const
  {
     return Anormal.norm();
  }

  // mean of vertex coordinates
  Vector Face::barycenter( void ) const
  {
     Vector c( 0., 0., 0. );
     double n = 0.;

     HalfEdgeCIter h = he;
     do
     {
        c += h->vertex->position;
        n += 1.;
        h = h->next;
     }
     while( h != he );

     return c/n;
  }

  // PS only needed by old code and dual code
  // return the 3-space vector corresponding to the distinguised direction
  const Vector Face::Xvector( void ) const {
    return X()->next->vertex->position - X()->vertex->position;
  }

  bool Face::isBoundary( void ) const { return he->onBoundary; }

  // which of the 3 possible choices is v? (used to figure out whether
  // a given vertex is i, j, or k (so to speak)
  unsigned int Face::i( const VertexCIter v ) const{
    const VertexCIter v0 = he->vertex, v1 = he->next->vertex, v2 = he->next->next->vertex;
    return v == v0 ? 0 : ( v == v1 ? 1 : ( assert( v == v2 ), 2 ) );
  }

  void Face::InitAreaAnglesAngleSums( void ){
    assert( !isBoundary() );
    assert( he->next->next->next == he ); // assuming triangles

    // masses
    A = Anormal.norm();

    // compute Euclidean angles and add to three corner slots
    const VertexIter vi = he->vertex, vj = he->next->vertex, vk = he->next->next->vertex;
    const Vector pi = vi->position, pj = vj->position, pk = vk->position;
    // normalized edges
    const Vector tij = pj-pi, tjk = pk-pj, tki = pi-pk;
    const double sinT = 2*A;
    vi->s += ( alpha[0] = atan2( sinT, -dot( tij, tki ) ) );
    vj->s += ( alpha[1] = atan2( sinT, -dot( tjk, tij ) ) );
    vk->s += ( alpha[2] = atan2( sinT, -dot( tki, tjk ) ) );
  }

  double Face::K( void ) const{
    const VertexCIter vi = he->vertex, vj = he->next->vertex, vk = he->next->next->vertex;
    return alpha[0] * ( vi->s - 1 ) + alpha[1] * ( vj->s - 1 ) + alpha[2] * ( vk->s - 1 );
  }

  Vector Face :: sampleUniform( void ) const
  {
     // [This subroutine dedicated to Jim Arvo.]
     
     // Compute the warping function (xi1,xi2) -> (s,t).
     double xi1 = unitRand();
     double xi2 = unitRand();
     double s = sqrt( xi1 );
     double t = xi2;

     // Plug the warped coords into the original parameterization.
     Vector A( 1., 0., 0. );
     Vector B( 0., 1., 0. );
     Vector C( 0., 0., 1. );
     Vector P = (1.-s)*A + s*(1.-t)*B + s*t*C;

     return P;
  }

  Vector Face::sampleField( const Vector& p, double n ) const
  {
     using namespace DDGConstants;

     // TODO handle singular faces

     // get incident halfedges
     const HalfEdgeCIter h[3] = {
        he,
        he->next,
        he->next->next
     };
     
     // get transport angles along each edge // REMOVE w
     Vector r( fmodPI(n*(h[0]->r()+h[0]->w())),
	       fmodPI(n*(h[1]->r()+h[1]->w())),
	       fmodPI(n*(h[2]->r()+h[2]->w())) );

     // compute associated total curvature for the triangle
     double Omega = r[0] + r[1] + r[2];

     // compute the right-hand side y
     Vector y = Omega*p - r;

     // compute the transport angles from each vertex to p
     Vector x( y[2] - y[0],
               y[0] - y[1],
               y[1] - y[2] );
     x /= 3.;

     // // get field vectors
     // Vector N = normal();
     // Quaternion Z[3];
     // for( int i = 0; i < 3; i++ )
     // {
     //    Z[i] = h[i]->vertex->kVec();
     //    double a = Z[i].norm();
     //    Z[i].im() -= dot(Z[i].im(),N)*N;
     //    Z[i].normalize();
     //    Z[i] *= a;
     // }

     // get field vectors
     Vector N = normal;
     Quaternion Z[3];
     for( int i = 0; i < 3; i++ )
     {
        HalfEdgeCIter X = h[i]->vertex->X();
        double s = h[i]->vertex->s;
        double theta = 0.;
        HalfEdgeCIter Y = h[i];
        do
        {
           Y = Y->flip->next;
           theta += Y->cornerAngle() * s;
        }
        while( Y->flip->next != X );

        double phi = h[i]->vertex->u.arg() - theta;
        double a = h[i]->vertex->u.norm();
        Quaternion w = ( h[i]->flip->vertex->position -
                         h[i]->vertex->position ).unit();
        Quaternion q( cos(phi/2.), sin(phi/2.)*N );
        Z[i] = a * (q*w*q.conj()).im();
     }

     // compute in-plane rotations by each transport angle
     Quaternion q[3];
     q[0] = Quaternion( cos(x[0]/2.), -sin(x[0]/2.)*N );
     q[1] = Quaternion( cos(x[1]/2.), -sin(x[1]/2.)*N );
     q[2] = Quaternion( cos(x[2]/2.), -sin(x[2]/2.)*N );


     // transport field vectors from vertices to p
     Quaternion Zp = p[0]*( q[0].conj()*Z[0]*q[0] ) +
                     p[1]*( q[1].conj()*Z[1]*q[1] ) +
                     p[2]*( q[2].conj()*Z[2]*q[2] ) ;

     return Zp.im();
  }

  Vector Face :: toWorldCoordinates( const Vector& b ) const
  {
     Vector pi = he->vertex->position;
     Vector pj = he->next->vertex->position;
     Vector pk = he->next->next->vertex->position;

     return b[0]*pi +
            b[1]*pj +
            b[2]*pk ;
  }
}

