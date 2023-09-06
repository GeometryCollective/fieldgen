#ifndef ALIASTABLE_H
#define ALIASTABLE_H

#include <vector>
#include <stack>

class AliasTable : public std::vector<double>
{
   public:
      // INTERFACE:
      // void build( void );
      // int sample( void ) const;

      // IMPLEMENTATION:
      void build( void )
      {
         const std::vector<double>& value( *this );
         std::stack<int> rich, poor;
         double mean = 0.;

         for( size_t i = 0; i < size(); i++ )
         {
            mean += value[i];
         }
         mean /= (double) size();

         alias.resize( size() );
         percent.resize( size() );

         for( size_t i = 0; i < size(); i++ )
         {
            alias[i] = i;
            percent[i] = 1e-8 + value[i] / mean;
            if( percent[i] < 1. )
            {
               poor.push( i );
            }
            else
            {
               rich.push( i );
            }
         }

         while( !poor.empty() )
         {
            int p = poor.top(); poor.pop();
            int r = rich.top();

            percent[r] -= ( 1. - percent[p] );
            alias[p] = r;

            if( percent[r] < 1. )
            {
               rich.pop();
               poor.push( r );
            }
         }
      }

      int sample( void ) const
      {
         int i = rand() % size();

         if( unitRand() < percent[i] )
         {
            return i;
         }

         return alias[i];
      }

   protected:
      double unitRand( void ) const
      {
         const double rRandMax = 1. / (double) RAND_MAX;

         return (double) rand() * rRandMax;
      }

      std::vector<int> alias;
      std::vector<double> percent;
};

#endif
