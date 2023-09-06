#include "Edge.h"
#include "Mesh.h"
#include "HalfEdge.h"

namespace DDG
{
  // only needed for old code
  const Vector Edge::Xvector( void ) const {
    return X()->next->vertex->position - X()->vertex->position;
  }

  // standard cotan Laplace coeff for this edge
  const double Edge::cot( void ) const {
    return .5 * ( ( he->onBoundary ? 0 : he->cot() ) +
		  ( he->flip->onBoundary ? 0 : he->flip->cot() ) );
  }
}

