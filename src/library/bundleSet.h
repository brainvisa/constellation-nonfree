#ifndef CONSTELLATION_BUNDLESET_H
#define CONSTELLATION_BUNDLESET_H

#include <connectomist/fibertracking/bundles.h>

namespace constel
{
  typedef std::vector< comist::FiberPoint > Fiber;
  typedef std::vector< Fiber > Fibers;
  typedef std::vector< Fibers > BundlesSet;
}

#endif

