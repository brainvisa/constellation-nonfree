#ifndef CONSTELLATION_SELECTBUNDLESFROMNAMES_H
#define CONSTELLATION_SELECTBUNDLESFROMNAMES_H

#include <aims/fibers/bundles.h>
#include <regex.h>


namespace constel {

  //--------------------------
  //  SelectBundlesFromNames
  //--------------------------
  class SelectBundlesFromNames
    : public aims::BundleProducer, public aims::BundleListener {

   public:
    typedef std::vector< aims::FiberPoint > Fiber;
    /** If use_fiber_names is not set, filtering is performed (as normal)
        using BundleInfo names, at bundleStarted() time.
        If use_fiber_names is set, filtering is performed at fiberTerminated()
        time. This is a bit unexpected, but name info is sometimes set when
        a fiber is entirely parsed (see SelectFiberListenerFromMesh)
    */
    SelectBundlesFromNames(std::vector< std::string > &select_bundles_name,
                           bool verbose = true, bool as_regex=false,
                           bool use_fiber_names = false);
    SelectBundlesFromNames();
    virtual ~SelectBundlesFromNames();

   protected:
    virtual void bundleStarted( const aims::BundleProducer &,
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

    std::set<std::string> _select_bundles_name;
    std::set<regex_t *> _regex;
    bool _bundle_selected;
    bool _verbose;
    bool _as_regex;
    bool _use_fiber_names;
    bool _new_bundle;
    std::string _last_name;
    int _last_id;
    Fiber _fiber;
  };

} // namespace constel

#endif // ifndef CONSTELLATION_SELECTBUNDLESFROMNAMES_H

