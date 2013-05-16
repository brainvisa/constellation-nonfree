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

    SelectBundlesFromNames( std::vector< std::string > &select_bundles_name,
                            bool verbose = true, bool as_regex=false );
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
  };

} // namespace constel

#endif // ifndef CONSTELLATION_SELECTBUNDLESFROMNAMES_H

