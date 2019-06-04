#ifndef DDG_SECTIONINTEGRALS_H
#define DDG_SECTIONINTEGRALS_H

#include "Complex.h"

extern DDG::Complex
DirichletIJ( const double s, const double gii, const double gij, const double gjj );

extern double
DirichletII( const double s, const double gjj, const double gjk, const double gkk );

extern DDG::Complex
MassIJ( const double s );

extern double MassII( void );

#endif /* SECTIONINTEGRALS */
