#ifndef DDG_UTILITY_H
#define DDG_UTILITY_H

#include <cstdlib>
#include <cmath>
#include <sys/time.h>
#include <execinfo.h>
#include <iostream>
#include "Complex.h"

namespace DDGConstants
{
   static DDG::Complex ii( 0., 1. );
   const double PI = 3.14159265358979323846;
}

namespace DDG
{
   inline double sqr( double x )
   {
      return x*x;
   }

   inline double unitRand( void )
   {
      const double rRandMax = 1. / (double) RAND_MAX;

      return rRandMax * (double) rand();
   }

   inline double seconds( int t0, int t1 )
   {
      return (double)(t1-t0) / (double) CLOCKS_PER_SEC;
   }

   inline const double fmodPI( const double a )
   {
      using namespace DDGConstants;
      return a - (2*PI) * floor( (a+PI) / (2*PI) );
   }

   inline double wallClock( void )
   {
     struct timeval t;
     struct timezone tz;
     gettimeofday( &t, &tz );
     return (double) t.tv_sec + 1e-6*(double) t.tv_usec;
   }

   inline void printTiming( const char* message, double t )
   {
      using namespace std;

      void* callstack[128];
      int stackDepth = backtrace( callstack, 128 );

      for( int i = 0; i < stackDepth; i++ )
      {
         cerr << " ";
      }
      cerr << "Wall clock time to " << message << ": " << t << " seconds." << endl;
   }
}

#endif
