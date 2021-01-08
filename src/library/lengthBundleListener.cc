#include <constellation/lengthBundleListener.h>

using namespace std;
using namespace carto;
using namespace aims;

namespace constel {

  LengthBundleListener::LengthBundleListener(string & fibersLengthFileName)
      : _fileName(fibersLengthFileName) {
    _file.open(_fileName.c_str(), fstream::out);
    if (_file.is_open()) cout << "FibersLength file opened: " << _fileName
      << endl;
  }

  LengthBundleListener::~LengthBundleListener() {
    _file.close();
  }

  void LengthBundleListener::bundleStarted(
      const BundleProducer &, const BundleInfo & /* bundleInfo */) {}

  void LengthBundleListener::fiberStarted(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    _fiberLength = 0;
    _fiberPointCount = 0;
  }

  void LengthBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint &fiberPoint) {
    ++_fiberPointCount;
    if (_fiberPointCount > 1) {
      _fiberLength += (fiberPoint - _antFiberPoint).norm();
    }
    _antFiberPoint = fiberPoint;
  }

  void LengthBundleListener::fiberTerminated(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    _file << int(_fiberLength) << endl;
  }
} // namespace constel

