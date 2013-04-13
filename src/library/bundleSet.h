#ifndef CONSTELLATION_BUNDLESET_H
#define CONSTELLATION_BUNDLESET_H

#include <connectomist/fibertracking/bundles.h>

namespace comist {
  typedef std::vector< comist::FiberPoint > Fiber;
  typedef std::vector< Fiber > Faisceaux;
  typedef std::vector< Faisceaux > BundlesSet;
}

#endif

