#ifndef CONSTELLATION_BUNDLESET_H
#define CONSTELLATION_BUNDLESET_H

#include <aims/fibers/bundles.h>

namespace constel {
  typedef std::vector<aims::FiberPoint> Fiber;
  typedef std::vector<Fiber> Fibers;
  typedef std::vector<Fibers> BundlesSet;
}

#endif

