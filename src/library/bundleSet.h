#ifndef CONSTELLATION_BUNDLESET_H
#define CONSTELLATION_BUNDLESET_H

#include <aims/fibers/bundles.h>

namespace constel {
  typedef aims::FiberPoint Point;
  typedef std::vector<Point> Fiber;
  typedef std::vector<Fiber> Fibers;
  typedef std::vector<Fibers> BundlesSet;


  class WeightedFiber
  {
    public:
      inline WeightedFiber();
      inline WeightedFiber( Fiber fiber );
      inline WeightedFiber( Fiber fiber, double weight);
      inline Fiber fiber() const;
      inline double weight() const;

    protected:
      Fiber _fiber;
      double _weight;
  };

  typedef std::vector<WeightedFiber> WeightedFibers;

  inline WeightedFiber::WeightedFiber()
  : _fiber( Fiber() ), _weight( 1.0 )
  {}

  inline WeightedFiber::WeightedFiber( Fiber fiber )
  : _fiber( fiber ), _weight( 1.0 )
  {}

  inline WeightedFiber::WeightedFiber( Fiber fiber, double weight )
  : _fiber( fiber ), _weight( weight )
  {}

  inline Fiber WeightedFiber::fiber() const
  {
    return _fiber;
  }

  inline double WeightedFiber::weight() const
  {
    return _weight;
  }

}

#endif
