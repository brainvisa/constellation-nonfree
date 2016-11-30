#include <constellation/bundleTools.h>
using namespace aims;
using namespace std;

namespace constel {

  FiberNameAccordingToMeshIntersectionBundleListener
    ::FiberNameAccordingToMeshIntersectionBundleListener(
        BundleInteractionReader &bundleInteractionReader,
        const string &fileName) : _fileName(fileName) {
    _bundleInteractionReader = &bundleInteractionReader;
    _file.open(_fileName.c_str(), fstream::out);
    if (_file.is_open())
      cout << "MeshIntersection storing file opened: " << _fileName << endl;
  }

  FiberNameAccordingToMeshIntersectionBundleListener
      ::~FiberNameAccordingToMeshIntersectionBundleListener() {
    _file.close();
  }

  void FiberNameAccordingToMeshIntersectionBundleListener::fiberTerminated(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    set<int> intersectionMeshIdentity_set;
    vector<int> fiberMeshIntersectionMeshIdentityVector
      = _bundleInteractionReader
        ->_listenedFiberInfo.getFiberMeshIntersectionMeshIdentityVector();

    for (uint meshId = 0;
         meshId < fiberMeshIntersectionMeshIdentityVector.size(); ++meshId) {
      intersectionMeshIdentity_set.insert(
          fiberMeshIntersectionMeshIdentityVector[meshId]);
    }
    if (not intersectionMeshIdentity_set.empty()) {
      set<int>::const_iterator set_it, set_itb, set_ite;
      set_itb = intersectionMeshIdentity_set.begin();
      set_ite = intersectionMeshIdentity_set.end(); 
      for (set_it = set_itb; set_it != set_ite; ++ set_it) {
        if (set_it == set_itb) {
          _file << *set_it << flush;
        } else {
          _file << "_" << *set_it << flush;
        }
      }
      _file << endl;
    } else {  //intersectionMeshIdentity_set is empty()
      _file << "trash" << endl;
    }
  }
} //namespace constel
