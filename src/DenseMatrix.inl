#include <cassert>
#include <cmath>

#include <algorithm>
#include <iostream>
using namespace std;

#include "DenseMatrix.h"
#include "LinearContext.h"
#include "Quaternion.h"
#include "SparseMatrix.h"
#include "Utility.h"

namespace DDG
{
   extern LinearContext context;

   template <class T>
   DenseMatrix<T> :: DenseMatrix( int m_, int n_ )
   // initialize an mxn matrix
   : m( m_ ),
     n( n_ ),
     cData( NULL )
   {
      data.resize( m*n );
      zero();
   }

   template <class T>
   DenseMatrix<T> :: DenseMatrix( const DenseMatrix<T>& A )
   // copy constructor
   : cData( NULL )
   {
      *this = A;
   }

   template <class T>
   DenseMatrix<T> :: ~DenseMatrix( void )
   // destructor
   {
      if( cData != NULL )
      {
         cholmod_l_free_dense( &cData, context );
      }
   }

   template <class T>
   DenseMatrix<T> DenseMatrix<T> :: transpose( void ) const
   {
      const DenseMatrix<T>& A( *this );
      DenseMatrix<T> AT( n, m );

      for( int i = 0; i < n; i++ )
      for( int j = 0; j < m; j++ )
      {
         AT(i,j) = A(j,i).conj();
      }

      return AT;
   }

   template <class T>
   SparseMatrix<T> DenseMatrix<T>::sparse( void )
   // converts to a sparse matrix
   {
      SparseMatrix<T> B;

      B = cholmod_l_dense_to_sparse( this->to_cholmod(), true, context );

      return B;
   }

   template <class T>
   DenseMatrix<T> DenseMatrix<T> :: operator*( const DenseMatrix<T>& B ) const
   // returns product of this matrix with B
   {
      const DenseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nColumns() == B.nRows() );

      DenseMatrix<T> AB( A.nRows(), B.nColumns() );

      for( int i = 0; i < A.nRows(); i++ )
      for( int j = 0; j < B.nColumns(); j++ )
      for( int k = 0; k < A.nColumns(); k++ )
      {
         AB( i, j ) += A( i, k ) * B( k, j );
      }

      return AB;
   }

   template <class T>
   void DenseMatrix<T> :: operator*=( const T& c )
   {
      DenseMatrix<T>& A( *this );

      for( int i = 0; i < m; i++ )
      for( int j = 0; j < n; j++ )
      {
         A(i,j) *= c;
      }
   }

   template <class T>
   void DenseMatrix<T> :: operator/=( const T& c )
   {
      DenseMatrix<T>& A( *this );

      for( int i = 0; i < m; i++ )
      for( int j = 0; j < n; j++ )
      {
         A(i,j) /= c;
      }
   }

   template <class T>
   DenseMatrix<T> DenseMatrix<T> :: operator+( const DenseMatrix<T>& B ) const
   // returns sum of this matrix with B
   {
      const DenseMatrix& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      DenseMatrix<T> C( nRows(), nColumns() );

      for( int i = 0; i < nRows(); i++ )
      for( int j = 0; j < nColumns(); j++ )
      {
         C(i,j) = A(i,j) + B(i,j);
      }

      return C;
   }

   template <class T>
   void DenseMatrix<T> :: operator+=( const DenseMatrix<T>& B )
   {
      DenseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      for( int i = 0; i < nRows(); i++ )
      for( int j = 0; j < nColumns(); j++ )
      {
         A(i,j) += B(i,j);
      }
   }

   template <class T>
   DenseMatrix<T> DenseMatrix<T> :: operator-( const DenseMatrix<T>& B ) const
   // returns difference of this matrix with B
   {
      const DenseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      DenseMatrix C( nRows(), nColumns() );

      for( int i = 0; i < nRows(); i++ )
      for( int j = 0; j < nColumns(); j++ )
      {
         C(i,j) = A(i,j) - B(i,j);
      }

      return C;
   }

   template <class T>
   void DenseMatrix<T> :: operator-=( const DenseMatrix<T>& B )
   {
      DenseMatrix<T>& A( *this );

      // make sure matrix dimensions agree
      assert( A.nRows() == B.nRows() );
      assert( A.nColumns() == B.nColumns() );

      for( int i = 0; i < nRows(); i++ )
      for( int j = 0; j < nColumns(); j++ )
      {
         A(i,j) -= B(i,j);
      }
   }

   template <class T>
   DenseMatrix<T> operator*( const T& c, const DenseMatrix<T>& A )
   {
      DenseMatrix<T> cA = A;

      cA *= c;

      return cA;
   }

   template <class T>
   DenseMatrix<T> operator*( const DenseMatrix<T>& A, double c )
   {
      return c*A;
   }

   template <class T>
   DenseMatrix<T> operator/( const DenseMatrix<T>& A, double c )
   {
      DenseMatrix<T> Ac = A;

      Ac /= c;

      return Ac;
   }

   template <class T>
   const DenseMatrix<T>& DenseMatrix<T> :: operator=( const DenseMatrix<T>& B )
   // copies B
   {
      if( cData )
      {
         cholmod_l_free_dense( &cData, context );
         cData = NULL;
      }

      m = B.m;
      n = B.n;
      data = B.data;

      return *this;
   }

   template <class T>
   int DenseMatrix<T> :: nRows( void ) const
   // returns the number of rows
   {
      return m;
   }

   template <class T>
   int DenseMatrix<T> :: nColumns( void ) const
   // returns the number of columns
   {
      return n;
   }

   template <class T>
   int DenseMatrix<T> :: length( void ) const
   // returns the size of the largest dimension
   {
      return max( m, n );
   }

   template <class T>
   void DenseMatrix<T> :: zero( const T& val )
   // sets all elements to val
   {
      for( int i = 0; i < m*n; i++ )
      {
         data[i] = val;
      }
   }

   template <class T>
   double DenseMatrix<T> :: norm( void ) const
   // returns the maximum magnitude of any entry
   {
      double maxNorm = 0.;

      for( int i = 0; i < m*n; i++ )
      {
         maxNorm = max( maxNorm, data[i].norm() );
      }

      return maxNorm;
   }

   template <class T>
   T& DenseMatrix<T> :: operator()( int row, int col )
   {
      assert( 0 <= row && row < m && 0 <= col && col < n );
      return data[row+m*col];
   }

   template <class T>
   T DenseMatrix<T> :: operator()( int row, int col ) const
   {
      assert( 0 <= row && row < m && 0 <= col && col < n );
      return data[row+m*col];
   }

   template <class T>
   T& DenseMatrix<T> :: operator()( int index )
   {
      assert( 0 <= index && index < m*n );
      return data[index];
   }

   template <class T>
   T DenseMatrix<T> :: operator()( int index ) const
   {
      assert( 0 <= index && index < m*n );
      return data[index];
   }

   template <class T>
   void DenseMatrix<T>::normalize( void )
   {
      double rFrobeniusNorm = 0.;

      for( int i = 0; i < m*n; i++ )
      {
         rFrobeniusNorm += data[i].norm2();
      }

      rFrobeniusNorm = 1./sqrt( rFrobeniusNorm );

      for( int i = 0; i < m*n; i++ )
      {
         data[i] *= rFrobeniusNorm;
      }
   }

   template <class T>
   T DenseMatrix<T>::sum( void ) const
   // returns the sum of all entries
   {
      T total( 0., 0. );

      for( int i = 0; i < m*n; i++ )
      {
         total += data[i];
      }

      return total;
   }

   template <class T>
   void DenseMatrix<T> :: removeMean( void )
   {
      T mean = 0.;
      int N = m*n;

      for( int i = 0; i < N; i++ )
      {
         mean += data[i];
      }

      mean /= (double) N;

      for( int i = 0; i < N; i++ )
      {
         data[i] -= mean;
      }
   }

   template <class T>
   double dot( const DenseMatrix<T>& x, const DenseMatrix<T>& y )
   // returns Euclidean inner product of x and y
   {
      return ( x.transpose() * y )(0);
   }

   template <class T>
   DenseMatrix<T> DenseMatrix<T>::operator-( void ) const
   // returns additive inverse of this matrix
   {
      const DenseMatrix<T>& A( *this );
      DenseMatrix<T> B( m, n );

      for( int i = 0; i < m; i++ )
      for( int j = 0; j < n; j++ )
      {
         B( i, j ) = -A( i, j );
      }

      return B;
   }
}

