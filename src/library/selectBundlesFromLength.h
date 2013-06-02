#ifndef CONSTELLATION_SELECTBUNDLESFROMLENGTH_H
#define CONSTELLATION_SELECTBUNDLESFROMLENGTH_H

#include <connectomist/fibertracking/bundles.h>


namespace constel
{

  //------------------//
  // SelectBundlesFromLength //
  //------------------//
  class SelectBundlesFromLength
    : public comist::BundleProducer, public comist::BundleListener
  {
  public:
    typedef std::vector< comist::FiberPoint > Fiber;

    SelectBundlesFromLength( float lmin, float lmax, bool verbose = true );
    SelectBundlesFromLength();
    virtual ~SelectBundlesFromLength();

  protected:


    virtual void bundleStarted( const comist::BundleProducer &,
                                const comist::BundleInfo &bundleInfo );
    virtual void bundleTerminated( const comist::BundleProducer &,
                                   const comist::BundleInfo &bundleInfo );
    virtual void fiberStarted( const comist::BundleProducer &,
                               const comist::BundleInfo &,
                               const comist::FiberInfo & );
    virtual void fiberTerminated( const comist::BundleProducer &,
                                  const comist::BundleInfo &bundleInfo,
                                  const comist::FiberInfo &fiberInfo );
    virtual void newFiberPoint( const comist::BundleProducer &,
                                const comist::BundleInfo &,
                                const comist::FiberInfo &,
                                const comist::FiberPoint & );
    virtual void noMoreBundle( const comist::BundleProducer & );

    float _lmin;
    float _lmax;
    bool _fiber_selected;
    bool _verbose;
    float _fiberLength;//in mm
    Fiber _fiber;
    comist::FiberPoint _antFiberPoint;
  };

} // namespace constel

#endif // ifndef CONSTELLATION_SELECTBUNDLESFROMLENGTH_H

