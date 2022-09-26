#ifndef CONSTELLATION_FIBERSWEIGHTSREADER_H
#define CONSTELLATION_FIBERSWEIGHTSREADER_H

#include <aims/fibers/bundles.h>


namespace constel {

  //---------------------------
  //  FibersWeightsReader
  //---------------------------
  class FibersWeightsReader
    : public aims::BundleProducer, public aims::BundleListener {

   public:
    FibersWeightsReader(std::string weightsFilename,
                        bool verbose = true);
    FibersWeightsReader();
    virtual ~FibersWeightsReader();

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
    std::string _weightsFilename;
    std::ifstream _weightsFile;
  };

} // namespace constel

#endif // ifndef CONSTELLATION_FIBERSWEIGHTSREADER_H
