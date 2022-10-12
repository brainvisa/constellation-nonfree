#ifndef CONSTELLATION_BUNDLESET_H
#define CONSTELLATION_BUNDLESET_H

#include <aims/fibers/bundles.h>

namespace constel {
  typedef aims::FiberPoint Point;
  typedef std::vector<Point> Fiber;
  typedef std::vector<Fiber> Fibers;
  typedef std::vector<Fibers> BundlesSet;
  typedef std::pair<Fiber, double> WeightedFiber;
  typedef std::vector<WeightedFiber> WeightedFibers;
}

#endif
