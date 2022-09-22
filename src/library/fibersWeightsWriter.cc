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
    string weightsFilename, bool verbose)
    : _weightsFilename(weightsFilename), _verbose(verbose) {}


  //---------------------------------------------------------------------------
  FibersWeightsWriter::~FibersWeightsWriter() {}


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::bundleStarted(
      const BundleProducer &, const BundleInfo &bundleInfo) {
    // string line;
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
      string s;
      vector<float> weights;
      getline(_weightsFile, s, ' ' );
      _weight = stof(s);
      weights.push_back( _weight );
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
    // fiberInfo.weight() = _weights[fiberInfo.id()];
    if ( _verbose ) {
      cout << "Fiber id: " << fiberInfo.id() << endl;
      // cout << "FiberInfo weight: " << fiberInfo.weight() <<endl;
      cout << "Fiber weight: " << _weight << endl;
    }
    // _fiber.clear();
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::fiberTerminated( const BundleProducer &,
                                     const BundleInfo &bundleInfo,
                                     const FiberInfo &fiberInfo ) {
      terminateFiber(bundleInfo, fiberInfo);
  }


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::newFiberPoint(
      const BundleProducer &, const BundleInfo & /* bundleInfo */,
      const FiberInfo & /* fiberInfo */, const FiberPoint &point) {}


  //---------------------------------------------------------------------------
  void FibersWeightsWriter::noMoreBundle(const BundleProducer &) {
      _weightsFile.close();
      BundleProducer::noMoreBundle();
      if (_verbose) cout << "end of FibersWeightsWriter filtering."
        << endl << flush;
  }


}
