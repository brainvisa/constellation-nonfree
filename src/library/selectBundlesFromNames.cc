#include <constellation/selectBundlesFromNames.h>

using namespace std;
using namespace carto;
using namespace aims;
using namespace comist;


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
SelectBundlesFromNames::SelectBundlesFromNames( vector< string > &select_bundles_name, bool verbose, bool as_regex ) :
    _verbose(verbose), _as_regex( as_regex )
{
  _select_bundles_name.insert( select_bundles_name.begin(),
                               select_bundles_name.end() );
  if (_verbose)
    for (set<string>::iterator i=_select_bundles_name.begin(),
         e=_select_bundles_name.end(); i!=e; i++)
      cout << *i << endl;

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


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::bundleStarted( const BundleProducer &,
                                   const BundleInfo &bundleInfo )
{
  if( _as_regex )
  {
    if( _regex.empty() )
    {
      // 1st time: compile regex for each string
      set<string>::iterator i, e = _select_bundles_name.end();
      for( i=_select_bundles_name.begin(); i!=e; ++i )
      {
        regex_t *preg = new regex_t;
        int status = regcomp( preg, i->c_str(), REG_EXTENDED );
        if( status == 0 )
          _regex.insert( preg );
        else
          delete preg;
      }
    }

    // try to match any regex to the bundle name
    regmatch_t pmatch;
    set<regex_t *>::iterator ir, er = _regex.end();
    int status = 1;
    for( ir=_regex.begin(); ir!=er && status!=0; ++ir )
      status &= regexec( *ir, bundleInfo.name().c_str(), 1, &pmatch, 0 );

    if( status == 0 )
      _bundle_selected = true;
    else
      _bundle_selected = false;
  }
  else
  {
    // recherche du nom du bundles
    set<string>::iterator recherche =
      _select_bundles_name.find(bundleInfo.name());
    _bundle_selected = !( recherche == _select_bundles_name.end() );
  }

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
  if ( _bundle_selected )
  {
    terminateBundle(bundleInfo);
  }
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::fiberStarted( const BundleProducer &,
                                  const BundleInfo &,
                                  const FiberInfo & )
{
  _fiber.clear();
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo )
{
  if ( _bundle_selected )
  {
    startFiber(bundleInfo, fiberInfo);
    for (size_t i = 0; i< _fiber.size(); i++)
    {
      addFiberPoint( bundleInfo, fiberInfo, _fiber[i]);
    }

    terminateFiber(bundleInfo, fiberInfo);
  }
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::newFiberPoint( const BundleProducer &, 
                                   const BundleInfo &,
                                   const FiberInfo &,
                                   const FiberPoint &point )
{
  _fiber.push_back(point);
}


//-----------------------------------------------------------------------------
void SelectBundlesFromNames::noMoreBundle( const BundleProducer & )
{
    BundleProducer::noMoreBundle();
    if (_verbose) cout << endl << "fin de l'écriture!!" << endl << flush;
}

} // namespace comist
