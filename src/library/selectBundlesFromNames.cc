#include <constellation/selectBundlesFromNames.h>

using namespace std;
using namespace carto;
using namespace aims;


namespace constel
{

//------------------//
// SelectBundlesFromNames //
//------------------//

//-----------------------------------------------------------------------------
SelectBundlesFromNames::SelectBundlesFromNames() :
    _verbose( false ), _as_regex( false )
{
}

//-----------------------------------------------------------------------------
SelectBundlesFromNames::SelectBundlesFromNames( vector< string > &select_bundles_name, bool verbose, bool as_regex, bool use_fiber_names ) :
    _verbose(verbose), _as_regex( as_regex ), 
    _use_fiber_names( use_fiber_names ), _new_bundle( true )
{
  _select_bundles_name.insert( select_bundles_name.begin(),
                               select_bundles_name.end() );
  if (_verbose)
    for (set<string>::iterator i=_select_bundles_name.begin(),
         e=_select_bundles_name.end(); i!=e; i++)
      cout << "selection: "<<*i << endl;
}


//-----------------------------------------------------------------------------
SelectBundlesFromNames::~SelectBundlesFromNames()
{
  set<regex_t *>::iterator i, e = _regex.end();
  for( i=_regex.begin(); i!=e; ++i )
  {
    regfree( *i );
    delete *i;
  }
}


namespace
{

  bool filterBundleName( const string & name, 
    const set<string> & select_bundles_name, bool as_regex, 
    set<regex_t *> & regex )
  {
    bool bundle_selected = false;

    if( as_regex )
    {
      if( regex.empty() )
      {
        // 1st time: compile regex for each string
        set<string>::iterator i, e = select_bundles_name.end();
        for( i=select_bundles_name.begin(); i!=e; ++i )
        {
          regex_t *preg = new regex_t;
          int status = regcomp( preg, i->c_str(), REG_EXTENDED );
          if( status == 0 )
            regex.insert( preg );
          else
            delete preg;
        }
      }

      // try to match any regex to the bundle name
      regmatch_t pmatch;
      set<regex_t *>::iterator ir, er = regex.end();
      int status = 1;
      for( ir=regex.begin(); ir!=er && status!=0; ++ir )
        status &= regexec( *ir, name.c_str(), 1, &pmatch, 0 );

      if( status == 0 )
        bundle_selected = true;
      else
        bundle_selected = false;
    }
    else
    {
      // search bundle name
      set<string>::iterator recherche =
        select_bundles_name.find( name );
      bundle_selected = !( recherche == select_bundles_name.end() );
    }

    return bundle_selected;
  }

}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::bundleStarted( const BundleProducer &,
                                            const BundleInfo &bundleInfo )
{
  if( _use_fiber_names )
  {
    // filter later, and also report startBundle()
    _bundle_selected = true;
    _new_bundle = true;
    return;
  }

  // here _use_fiber_names is not set: filter now
  _bundle_selected = filterBundleName( bundleInfo.name(), _select_bundles_name, 
                                       _as_regex, _regex );

  if( !_bundle_selected )
  {
    if (_verbose) cout << "bundle name: " << bundleInfo.name() << "... not selected" << endl;
  }
  else
  {
    startBundle(bundleInfo);
    if (_verbose) cout << "bundle name: " << bundleInfo.name() << "... OK" << endl;
  }
}

//-----------------------------------------------------------------------------
void SelectBundlesFromNames::bundleTerminated( const BundleProducer &,
                                      const BundleInfo &bundleInfo )
{
  if( _use_fiber_names && _bundle_selected )
  {
    BundleInfo binfo2( _last_id, _last_name );
    terminateBundle( binfo2 );
    return;
  }
  if ( _bundle_selected )
    terminateBundle(bundleInfo);
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::fiberStarted( const BundleProducer &,
                                  const BundleInfo & bundleInfo,
                                  const FiberInfo & fiberInfo )
{
  if ( _bundle_selected )
  {
    if( _use_fiber_names )
    {
      _fiber.clear();
    }
    else
      startFiber(bundleInfo, fiberInfo);
  }
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo )
{
  if ( _bundle_selected )
  {
    if( _use_fiber_names )
    {
      // select now
      bool selected = filterBundleName( bundleInfo.name(),
        _select_bundles_name, _as_regex, _regex );
      if( selected )
      {
        if( _new_bundle )
        {
          startBundle( bundleInfo );
          _new_bundle = false;
          _last_name = bundleInfo.name();
          _last_id = bundleInfo.id();
        }
        else if( bundleInfo.name() != _last_name )
        {
          // change/split bundle
          BundleInfo binfo( _last_id, _last_name );
          terminateBundle( binfo );
          startBundle( bundleInfo );
          _last_name = bundleInfo.name();
          _last_id = bundleInfo.id();
        }

        startFiber( bundleInfo, fiberInfo );
        Fiber::const_iterator ip, ep = _fiber.end();
        for( ip=_fiber.begin(); ip!=ep; ++ip )
          addFiberPoint( bundleInfo, fiberInfo, *ip );
        terminateFiber(bundleInfo, fiberInfo);
      }
    }
    else
      terminateFiber(bundleInfo, fiberInfo);
    _fiber.clear();
  }
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::newFiberPoint( const BundleProducer &, 
                                   const BundleInfo & bundleInfo,
                                   const FiberInfo & fiberInfo,
                                   const FiberPoint &point )
{
  if ( _bundle_selected )
  {
    if( _use_fiber_names )
      _fiber.push_back( point );
    else
      addFiberPoint( bundleInfo, fiberInfo, point );
  }
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::noMoreBundle( const BundleProducer & )
{
    BundleProducer::noMoreBundle();
    if (_verbose) cout << endl << "end of SelectBundlesFromNames filtering."
      << endl << flush;
}

} // namespace constel
