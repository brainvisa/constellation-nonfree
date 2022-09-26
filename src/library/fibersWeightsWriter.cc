#include <constellation/fibersWeightsWriter.h>

using namespace std;
using namespace carto;
using namespace aims;

namespace constel {

  //---------------------------------------------------------------------------
  FibersWeightsWriter::FibersWeightsWriter() :
    _verbose(false) {}

  //---------------------------------------------------------------------------
  FibersWeightsWriter::FibersWeightsWriter(
    string outputWeightsFilename,  bool verbose)
    : _outputWeightsFilename(outputWeightsFilename),
      _verbose(verbose) {}


  //---------------------------------------------------------------------------
  FibersWeightsWriter::~FibersWeightsWriter() {}


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::bundleStarted(
      const BundleProducer &, const BundleInfo &bundleInfo) {

    // Open output weights file
    _outputWeightsFile.open( _outputWeightsFilename );
    if ( _outputWeightsFile.fail() ) {
      if ( _verbose ) {
        cout << "Problem opening file: " << _outputWeightsFilename << endl;
      }
      exit(1);
    }
    else {
      if ( _verbose ) {
        cout << "Opened file: " << _outputWeightsFilename << endl;
      }
    }

    startBundle(bundleInfo); // needed?
  }

  //---------------------------------------------------------------------------
  void FibersWeightsWriter::bundleTerminated(
      const BundleProducer &, const BundleInfo &bundleInfo) {
    if (bundleInfo.name() == "255")
      cout << "FibersWeightsWriter, bundle 255 found, id: "
      << bundleInfo.id() << endl;
    terminateBundle(bundleInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::fiberStarted(
      const BundleProducer &, const BundleInfo & /* bundleInfo */,
      const FiberInfo &fiberInfo /* fiberInfo */) {
    if ( _verbose ) {
      cout << "Begin fiber id: " << fiberInfo.id() << endl;
    }

    // _fiber.clear();
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo ) {
      float weight = fiberInfo.weight();
      _outputWeightsFile << weight << ' ';

      if ( _verbose ) {
        cout << "[FibersWeightsWriter] Ending fiber terminated (id, weight): "
        << fiberInfo.id() << " , " << weight << endl;
      }

      terminateFiber(bundleInfo, fiberInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::newFiberPoint(
      const BundleProducer &, const BundleInfo & /* bundleInfo */,
      const FiberInfo & /* fiberInfo */, const FiberPoint &point) {}


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::noMoreBundle(const BundleProducer &) {
      _outputWeightsFile.close();
      BundleProducer::noMoreBundle();
      if (_verbose) cout << "end of FibersWeightsWriter filtering."
        << endl << flush;
  }


}
