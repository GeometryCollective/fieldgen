#include <cmath>
#include <iostream>
#include <assert.h>

using namespace std;

#include "Quaternion.h"
#include "Complex.h"

namespace DDG
{
   // CONSTRUCTORS ----------------------------------------------------------
   
   Quaternion :: Quaternion( void )
   // initializes all components to zero
   : s( 0. ),
     v( 0., 0., 0. )
   {}
   
   Quaternion :: Quaternion( const Quaternion& q )
   // initializes from existing quaternion
   : s( q.s ),
     v( q.v )
   {}
   
   Quaternion :: Quaternion( double s_, double vi, double vj, double vk )
   // initializes with specified double (s) and imaginary (v) components
   : s( s_ ),
     v( vi, vj, vk )
   {}
   
   Quaternion :: Quaternion( double s_, const Vector& v_ )
   // initializes with specified double(s) and imaginary (v) components
   : s( s_ ),
     v( v_ )
   {}
   
   Quaternion :: Quaternion( double s_ )
     : s( s_ ),
       v( 0., 0., 0. )
   {}
   
   Quaternion :: Quaternion( const Vector& v_ )
     : s( 0. ),
       v( v_ )
   {}

  Quaternion::Quaternion( const Complex& z )
    : s( z.re ), v( z.im, 0., 0. )
  {
  }
   
  // these codes are from SIGGRAPH course notes 23, Math for SIGGRAPH,
  // 1989, page 203, 204
  Quaternion::Quaternion( const Vector& e1, const Vector& e2, const Vector& e3 )
    : s(0.), v(0.,0.,0.)
  {
    assert( abs( e1.norm2() - 1 ) < 1e-12 );
    assert( abs( e2.norm2() - 1 ) < 1e-12 );
    assert( abs( e3.norm2() - 1 ) < 1e-12 );
    assert( dot( e1, e2 ) < 1e-12 );
    assert( dot( e2, e3 ) < 1e-12 );
    assert( dot( e3, e1 ) < 1e-12 );
    assert( abs( dot( cross( e1, e2 ), e3 ) - 1. ) < 1e-12 );

    // numerically stable way to convert a 3x3 matrix given by 3
    // orthonormal columns into a quaternion which affects that
    // rotation
    const double m[3][3] =
      { { e1[0], e2[0], e3[0] },
	{ e1[1], e2[1], e3[1] },
	{ e1[2], e2[2], e3[2] } };
    const double trace = m[0][0] + m[1][1] + m[2][2];
    if( trace > 0. ){
      double w = sqrt( trace + 1 );
      s = w * .5;
      w = .5 / w;
      v[0] = ( m[2][1] - m[1][2] ) * w;
      v[1] = ( m[0][2] - m[2][0] ) * w;
      v[2] = ( m[1][0] - m[0][1] ) * w;
    }else{
      unsigned int i = 0;
      if( m[1][1] > m[0][0] )	i = 1;
      if( m[2][2] > m[i][i] )	i = 2;
      const unsigned int j = ( i + 1 ) % 3;
      const unsigned int k = ( j + 1 ) % 3;
      double w = sqrt( ( m[i][i] - ( m[j][j] + m[k][k] ) ) + 1.0 );
      v[i] = w * .5;
      w = .5 / w;
      s = ( m[k][j] - m[j][k] ) * w;
      v[j] = ( m[j][i] + m[i][j] ) * w;
      v[k] = ( m[k][i] + m[i][k] ) * w;
    }
    // if( s < 0 ){ (*this) *= -1; cerr << "q"; }
  }
   
   // ASSIGNMENT OPERATORS --------------------------------------------------
   
   const Quaternion& Quaternion :: operator=( double _s )
   // assigns a purely real quaternion with real value s
   {
      s = _s;
      v = Vector( 0., 0., 0. );
   
      return *this;
   }
   
   const Quaternion& Quaternion :: operator=( const Vector& _v )
   // assigns a purely real quaternion with imaginary value v
   {
      s = 0.;
      v = _v;
   
      return *this;
   }
   
   
   // ACCESSORS -------------------------------------------------------------
   
   double& Quaternion::operator[]( int index )
   // returns reference to the specified component (0-based indexing: double, i, j, k)
   {
      return ( &s )[ index ];
   }
   
   const double& Quaternion::operator[]( int index ) const
   // returns const reference to the specified component (0-based indexing: double, i, j, k)
   {
      return ( &s )[ index ];
   }
   
   void Quaternion::toMatrix( double Q[4][4] ) const
   // returns 4x4 matrix representation
   {
      Q[0][0] =   s; Q[0][1] = -v.x; Q[0][2] = -v.y; Q[0][3] = -v.z;
      Q[1][0] = v.x; Q[1][1] =    s; Q[1][2] = -v.z; Q[1][3] =  v.y;
      Q[2][0] = v.y; Q[2][1] =  v.z; Q[2][2] =    s; Q[2][3] = -v.x;
      Q[3][0] = v.z; Q[3][1] = -v.y; Q[3][2] =  v.x; Q[3][3] =    s;
   }
   
   double& Quaternion::re( void )
   // returns reference to double part
   {
      return s;
   }
   
   const double& Quaternion::re( void ) const
   // returns const reference to double part
   {
      return s;
   }
   
   Vector& Quaternion::im( void )
   // returns reference to imaginary part
   {
      return v;
   }
   
   const Vector& Quaternion::im( void ) const
   // returns const reference to imaginary part
   {
      return v;
   }
   
   
   // VECTOR SPACE OPERATIONS -----------------------------------------------
   
   Quaternion Quaternion::operator+( const Quaternion& q ) const
   // addition
   {
      return Quaternion( s+q.s, v+q.v );
   }
   
   Quaternion Quaternion::operator-( const Quaternion& q ) const
   // subtraction
   {
      return Quaternion( s-q.s, v-q.v );
   }
   
   Quaternion Quaternion::operator-( void ) const
   // negation
   {
      return Quaternion( -s, -v );
   }
   
   Quaternion Quaternion::operator*( double c ) const
   // scalar multiplication
   {
      return Quaternion( s*c, v*c );
   }
   
   Quaternion operator*( double c, const Quaternion& q )
   // scalar multiplication
   {
      return q*c;
   }
   
   Quaternion Quaternion::operator/( double c ) const
   // scalar division
   {
      return Quaternion( s/c, v/c );
   }
   
   void Quaternion::operator+=( const Quaternion& q )
   // addition / assignment
   {
      s += q.s;
      v += q.v;
   }
   
   void Quaternion::operator+=( double c )
   // addition / assignment of pure real
   {
      s += c;
   }
   
   void Quaternion::operator-=( const Quaternion& q )
   // subtraction / assignment
   {
      s -= q.s;
      v -= q.v;
   }
   
   void Quaternion::operator-=( double c )
   // subtraction / assignment of pure real
   {
      s -= c;
   }
   
   void Quaternion::operator*=( double c )
   // scalar multiplication / assignment
   {
      s *= c;
      v *= c;
   }
   
   void Quaternion::operator/=( double c )
   // scalar division / assignment
   {
      s /= c;
      v /= c;
   }
   
   
   // ALGEBRAIC OPERATIONS --------------------------------------------------
   
   Quaternion Quaternion::operator*( const Quaternion& q ) const
   // Hamilton product
   {
      const double& s1( s );
      const double& s2( q.s );
      const Vector& v1( v );
      const Vector& v2( q.v );
   
      return Quaternion( s1*s2 - dot(v1,v2), s1*v2 + s2*v1 + cross(v1,v2) );
   }
   
   void Quaternion::operator*=( const Quaternion& q )
   // Hamilton product / assignment
   {
      *this = ( *this * q );
   }
   
   Quaternion Quaternion::conj( void ) const
   // conjugation
   {
      return Quaternion( s, -v );
   }
   
   Quaternion Quaternion::inv( void ) const
   {
      return ( this->conj() ) / this->norm2();
   }
   
   
   // NORMS -----------------------------------------------------------------
   
   double Quaternion::norm( void ) const
   // returns Euclidean length
   {
      return sqrt( s*s + v.x*v.x + v.y*v.y + v.z*v.z );
   }
   
   double Quaternion::norm2( void ) const
   // returns Euclidean length squared
   {
      return s*s + dot(v,v);
   }
   
   Quaternion Quaternion::unit( void ) const
   // returns unit quaternion
   {
      return *this / norm();
   }
   
   void Quaternion::normalize( void )
   // divides by Euclidean length
   {
      *this /= norm();
   }
   
   
   // GEOMETRIC OPERATIONS --------------------------------------------------
   
   Quaternion slerp( const Quaternion& q0, const Quaternion& q1, double t )
   // spherical-linear interpolation
   {
      // interpolate length
      double m0 = q0.norm();
      double m1 = q1.norm();
      double m = (1-t)*m0 + t*m1;
   
      // interpolate direction
      Quaternion p0 = q0 / m0;
      Quaternion p1 = q1 / m1;
      double theta = acos(( p0.conj()*p1 ).re() );
      Quaternion p = ( sin((1-t)*theta)*p0 + sin(t*theta)*p1 )/sin(theta);
   
      return m*p;
   }
   
   
   // I/O -------------------------------------------------------------------------
   
   std::ostream& operator<<( std::ostream& os, const Quaternion& q )
   // prints components
   {
      os << "( " << q.re() << ", " << q.im() << " )";
   
      return os;
   }
}

