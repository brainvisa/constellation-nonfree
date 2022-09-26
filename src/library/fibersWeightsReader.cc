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
      _verbose(verbose) {}


  //---------------------------------------------------------------------------
  FibersWeightsReader::~FibersWeightsReader() {}


  //---------------------------------------------------------------------------
  void FibersWeightsReader::bundleStarted(
      const BundleProducer &, const BundleInfo &bundleInfo) {
    // Open weights file
    _weightsFile.open( _weightsFilename );
    if ( _weightsFile.fail() ) {
      if ( _verbose ) {
        cout << "Problem reading file: " << _weightsFilename << endl;
      }
      exit(1);
    }
    else {
      if ( _verbose ) {
        cout << "Opened file: " << _weightsFilename << endl;
      }
    }

    startBundle(bundleInfo); // needed?
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
      const BundleProducer &, const BundleInfo & /* bundleInfo */,
      const FiberInfo &fiberInfo /* fiberInfo */) {
    if ( _verbose ) {
      cout << "Begin fiber id: " << fiberInfo.id() << endl;
    }

    // _fiber.clear();
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo ) {
      string weight;
      getline(_weightsFile, weight, ' ' );
      FiberInfo newFiberInfo( fiberInfo.id(), stof(weight) );

      if ( _verbose ) {
        cout << "[FibersWeightsReader] Ending fiber terminated (id, weight): "
        << fiberInfo.id() << " , " << weight << endl;
      }

      terminateFiber(bundleInfo, newFiberInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsReader::newFiberPoint(
      const BundleProducer &, const BundleInfo & /* bundleInfo */,
      const FiberInfo & /* fiberInfo */, const FiberPoint &point) {}


  //---------------------------------------------------------------------------
  void FibersWeightsReader::noMoreBundle(const BundleProducer &) {
      _weightsFile.close();
      BundleProducer::noMoreBundle();
      if (_verbose) cout << "end of FibersWeightsReader filtering."
        << endl << flush;
  }


}
