#include <constellation/fibersWeightsReader.h>

using namespace std;
using namespace carto;
using namespace aims;

namespace constel {

  //---------------------------------------------------------------------------
  FibersWeightsReader::FibersWeightsReader() :
    _verbose(false) {}

  //---------------------------------------------------------------------------
  FibersWeightsReader::FibersWeightsReader(
    string weightsFilename,  bool verbose)
    : _weightsFilename(weightsFilename),
      _verbose(verbose)
      {
        // Open weights file
        _weightsFile.open( _weightsFilename );
        if ( _weightsFile.fail() )
        {
          if ( _verbose )
          {
            cout << "Problem reading file: " << _weightsFilename << endl;
          }
          exit(1);
        }
        else
        {
          if ( _verbose )
          {
            cout << "Opened file: " << _weightsFilename << endl;
          }
        }
      }


  //---------------------------------------------------------------------------
  FibersWeightsReader::~FibersWeightsReader() {
    _weightsFile.close();
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::bundleStarted(
      const BundleProducer &, const BundleInfo &bundleInfo) {


    startBundle(bundleInfo);
  }

  //---------------------------------------------------------------------------
  void FibersWeightsReader::bundleTerminated(
      const BundleProducer &, const BundleInfo &bundleInfo) {
    if (bundleInfo.name() == "255")
      cout << "FibersWeightsReader, bundle 255 found, id: "
      << bundleInfo.id() << endl;
    terminateBundle(bundleInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::fiberStarted(
      const BundleProducer &, const BundleInfo &bundleInfo /* bundleInfo */,
      const FiberInfo &fiberInfo /* fiberInfo */) {
    string weight;
    getline(_weightsFile, weight, ' ' );
    if (weight == "#") {
      /* command history in weight file */
      getline(_weightsFile, weight, '\n');
      getline(_weightsFile, weight, ' ');

    }
    FiberInfo newFiberInfo( fiberInfo.id(), stod(weight) );
    _newFiberInfo = newFiberInfo;
    startFiber(bundleInfo, _newFiberInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo ) {
      terminateFiber(bundleInfo, _newFiberInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::newFiberPoint(
      const BundleProducer &, const BundleInfo &  bundleInfo,
      const FiberInfo & /* fiberInfo */, const FiberPoint &point) {
        addFiberPoint(bundleInfo, _newFiberInfo, point);
      }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::noMoreBundle(const BundleProducer &) {
      _weightsFile.close();
      BundleProducer::noMoreBundle();
      if (_verbose) cout << "end of FibersWeightsReader filtering."
        << endl << flush;
  }


}
