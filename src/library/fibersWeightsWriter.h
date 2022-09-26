#ifndef CONSTELLATION_FIBERSWEIGHTSWRITER_H
#define CONSTELLATION_FIBERSWEIGHTSWRITER_H

#include <aims/fibers/bundles.h>


namespace constel {

  //---------------------------
  //  FibersWeightsWriter
  //---------------------------
  class FibersWeightsWriter
    : public aims::BundleProducer, public aims::BundleListener {

   public:
    FibersWeightsWriter(std::string outputWeightsFilename,
                        bool verbose = true);
    FibersWeightsWriter();
    virtual ~FibersWeightsWriter();

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

  private:
    bool _verbose;
    std::string _outputWeightsFilename;
    std::ofstream _outputWeightsFile;
  };

} // namespace constel

#endif // ifndef CONSTELLATION_FIBERSWEIGHTSLISTENER_H
