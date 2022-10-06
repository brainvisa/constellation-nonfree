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
      _verbose(verbose) {
        // Open output weights file
        _outputWeightsFile.open( _outputWeightsFilename );
        if ( _outputWeightsFile.fail() ) {
          if ( _verbose ) {
            cout << "Problem opening file: " << _outputWeightsFilename << endl;
            cout << _outputWeightsFile.fail() << endl;
          }
          exit(1);
        }
        else {
          if ( _verbose ) {
            cout << "Opened file: " << _outputWeightsFilename << endl;
          }
        }
      }


  //---------------------------------------------------------------------------
  FibersWeightsWriter::~FibersWeightsWriter() {
    _outputWeightsFile << flush;
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::bundleStarted(
      const BundleProducer &, const BundleInfo &bundleInfo) {}

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
      const FiberInfo & /* fiberInfo */) {}


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo ) {
      _outputWeightsFile << std::setprecision(10) << fiberInfo.weight() << ' ';
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
