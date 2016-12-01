#ifndef CONSTELLATION_LENGTHBUNDLELISTENER_H
#define CONSTELLATION_LENGTHBUNDLELISTENER_H

#include <aims/fibers/bundles.h>

namespace constel {

  class LengthBundleListener : public aims::BundleListener {
   public:

    LengthBundleListener(std::string & fibersLengthFileName);

    virtual void bundleStarted(const aims::BundleProducer &,
                                const aims::BundleInfo &);
    virtual void bundleTerminated(const aims::BundleProducer &,
                                   const aims::BundleInfo &);
    virtual void fiberStarted(const aims::BundleProducer &,
                               const aims::BundleInfo &,
                               const aims::FiberInfo &);
    virtual void fiberTerminated(const aims::BundleProducer &,
                                  const aims::BundleInfo &,
                                  const aims::FiberInfo &);
    virtual void newFiberPoint(const aims::BundleProducer &,
                                const aims::BundleInfo &,
                                const aims::FiberInfo &,
                                const aims::FiberPoint &);
    virtual void noMoreBundle(const aims::BundleProducer &);

    virtual ~LengthBundleListener();

   protected:
    float _fiberLength; //in mm
    aims::FiberPoint _antFiberPoint;
    uint _fiberPointCount;
    std::string _fileName;
    std::fstream _file;
  };

  inline void LengthBundleListener::bundleTerminated(
      const aims::BundleProducer &, const aims::BundleInfo &) {}

  inline void LengthBundleListener::noMoreBundle(
      const aims::BundleProducer &) {}

} // namespace constel

#endif // ifndef CONSTELLATION_LENGTHBUNDLELISTENER_H
