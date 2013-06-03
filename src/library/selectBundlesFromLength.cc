#include <constellation/selectBundlesFromLength.h>

using namespace std;
using namespace carto;
using namespace aims;
using namespace comist;


namespace constel
{

//------------------//
// SelectBundlesFromLength //
//------------------//

//-----------------------------------------------------------------------------
SelectBundlesFromLength::SelectBundlesFromLength() :
    _lmin( 0 ), _lmax( -1 ), _verbose( false ), _fiberLength( 0 )
{
}

//-----------------------------------------------------------------------------
SelectBundlesFromLength::SelectBundlesFromLength( float lmin, float lmax,
  bool verbose ) :
    _lmin( lmin ), _lmax( lmax ), _verbose(verbose), _fiberLength( 0 )
{
  if (_verbose)
  {
    cout << "length min: " << _lmin << endl;
    cout << "length max: " << _lmax << endl;
  }
}


//-----------------------------------------------------------------------------
SelectBundlesFromLength::~SelectBundlesFromLength()
{
}


//-----------------------------------------------------------------------------
void SelectBundlesFromLength::bundleStarted( const BundleProducer &,
                                   const BundleInfo &bundleInfo )
{
  startBundle( bundleInfo );
}

//-----------------------------------------------------------------------------
void SelectBundlesFromLength::bundleTerminated( const BundleProducer &,
                                      const BundleInfo &bundleInfo )
{
  terminateBundle(bundleInfo);
}


//-----------------------------------------------------------------------------
void SelectBundlesFromLength::fiberStarted( const BundleProducer &,
                                  const BundleInfo & bundleInfo,
                                  const FiberInfo & fiberInfo )
{
  _fiberLength = 0;
  _fiber.clear();
}


//-----------------------------------------------------------------------------
void SelectBundlesFromLength::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo )
{
  if ( _fiberLength >= _lmin && ( _lmax < 0 || _fiberLength <= _lmax ) )
  {
//     if( _verbose )
//       cout << "fiber selected, length: " << _fiberLength << endl;
    startFiber( bundleInfo, fiberInfo );
    vector<FiberPoint>::const_iterator ip, ep = _fiber.end();
    for( ip=_fiber.begin(); ip!=ep; ++ip )
      addFiberPoint( bundleInfo, fiberInfo, *ip );
    terminateFiber( bundleInfo, fiberInfo );
  }
//   else if( _verbose )
//     cout << "fiber rejected, length: " << _fiberLength << endl;
  _fiber.clear();
  _fiberLength = 0;
}


//-----------------------------------------------------------------------------
void SelectBundlesFromLength::newFiberPoint( const BundleProducer &, 
                                   const BundleInfo & bundleInfo,
                                   const FiberInfo & fiberInfo,
                                   const FiberPoint &point )
{
  if( !_fiber.empty() )
    _fiberLength += ( point - _antFiberPoint ).norm();
  _antFiberPoint = point;
  _fiber.push_back( point );
}


//-----------------------------------------------------------------------------
void SelectBundlesFromLength::noMoreBundle( const BundleProducer & )
{
    BundleProducer::noMoreBundle();
    if (_verbose) cout << endl << "end of SelectBundlesFromLength filtering." 
      << endl << flush;
}

} // namespace comist
