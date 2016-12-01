#ifndef CONSTELLATION_SELECTBUNDLESFROMLENGTH_H
#define CONSTELLATION_SELECTBUNDLESFROMLENGTH_H

#include <aims/fibers/bundles.h>


namespace constel {

  //---------------------------
  //  SelectBundlesFromLength 
  //---------------------------
  class SelectBundlesFromLength
    : public aims::BundleProducer, public aims::BundleListener {

   public:
    typedef std::vector<aims::FiberPoint> Fiber;

    SelectBundlesFromLength(float lmin, float lmax, bool verbose = true);
    SelectBundlesFromLength();
    virtual ~SelectBundlesFromLength();

   protected:
    virtual void bundleStarted(const aims::BundleProducer &,
                                const aims::BundleInfo &bundleInfo);
    virtual void bundleTerminated(const aims::BundleProducer &,
                                   const aims::BundleInfo &bundleInfo);
    virtual void fiberStarted(const aims::BundleProducer &,
                               const aims::BundleInfo &,
                               const aims::FiberInfo &);
    virtual void fiberTerminated(const aims::BundleProducer &,
                                  const aims::BundleInfo &bundleInfo,
                                  const aims::FiberInfo &fiberInfo);
    virtual void newFiberPoint(const aims::BundleProducer &,
                                const aims::BundleInfo &,
                                const aims::FiberInfo &,
                                const aims::FiberPoint &);
    virtual void noMoreBundle(const aims::BundleProducer &);

    float _lmin;
    float _lmax;
    bool _fiber_selected;
    bool _verbose;
    float _fiberLength; //in mm
    Fiber _fiber;
    aims::FiberPoint _antFiberPoint;
  };

} // namespace constel

#endif // ifndef CONSTELLATION_SELECTBUNDLESFROMLENGTH_H

