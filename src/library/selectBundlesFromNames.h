#ifndef CONSTELLATION_SELECTBUNDLESFROMNAMES_H
#define CONSTELLATION_SELECTBUNDLESFROMNAMES_H

#include <connectomist/fibertracking/bundles.h>
#include <regex.h>


namespace constel
{

  //------------------//
  // SelectBundlesFromNames //
  //------------------//
  class SelectBundlesFromNames
    : public comist::BundleProducer, public comist::BundleListener
  {
  public:
    typedef std::vector< comist::FiberPoint > Fiber;

    /** If use_fiber_names is not set, filtering is performed (as normal)
        using BundleInfo names, at bundleStarted() time.
        If use_fiber_names is set, filtering is performed at fiberTerminated()
        time. This is a bit unexpected, but name info is sometimes set when
        a fiber is entirely parsed (see SelectFiberListenerFromMesh)
    */
    SelectBundlesFromNames( std::vector< std::string > &select_bundles_name,
                            bool verbose = true, bool as_regex=false,
                            bool use_fiber_names = false );
    SelectBundlesFromNames();
    virtual ~SelectBundlesFromNames();

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

