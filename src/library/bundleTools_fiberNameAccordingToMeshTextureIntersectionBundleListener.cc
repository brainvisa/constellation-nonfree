#include <constellation/bundleTools.h>

using namespace aims;
using namespace std;

namespace constel {

  FiberNameAccordingToMeshTextureIntersectionBundleListener
    ::FiberNameAccordingToMeshTextureIntersectionBundleListener(
        BundleInteractionReader &bundleInteractionReader,
        const TimeTexture<short> &labeled_tex,
        const string &fileName)
      : _fileName(fileName), _labeledTex0(labeled_tex.begin()->second) {
    _bundleInteractionReader = &bundleInteractionReader;
    _file.open(_fileName.c_str(), fstream::out);
    if (_file.is_open())
      cout << "MeshIntersection storing file opened: " << _fileName << endl;
  }

  FiberNameAccordingToMeshTextureIntersectionBundleListener
      ::~FiberNameAccordingToMeshTextureIntersectionBundleListener() {
    _file.close();
  }

  void FiberNameAccordingToMeshTextureIntersectionBundleListener::fiberStarted(
      const BundleProducer & /* bundleProducer */,
      const BundleInfo & /* bundleInfo */,
      const FiberInfo & /* fiberInfo */) {
    _intersectionPointsMeshClosestPointsTexLabels_set.clear();
  }

  void FiberNameAccordingToMeshTextureIntersectionBundleListener
      ::newFiberPoint(const BundleProducer & /* bundleProducer */,
                      const BundleInfo & /* bundleInfo */,
                      const FiberInfo & /* fiberInfo */,
                      const FiberPoint & /* point */ ) {
    if ((_bundleInteractionReader
         ->_listenedFiberInfo.getAntFiberPointExistingMeshIntersection())) {
      size_t antFiberPointMeshClosestPoint_index
        = _bundleInteractionReader
          ->_listenedFiberInfo.getAntFiberPointMeshClosestPointIndex();
      _intersectionPointsMeshClosestPointsTexLabels_set.insert(
          _labeledTex0[antFiberPointMeshClosestPoint_index]);
    }
  }


  void FiberNameAccordingToMeshTextureIntersectionBundleListener
      ::fiberTerminated(const BundleProducer &,
                        const BundleInfo &,
                        const FiberInfo &) {

    if (_bundleInteractionReader
        ->_listenedFiberInfo.getAntFiberPointExistingMeshIntersection()) {
      size_t antFiberPointMeshClosestPoint_index
        = _bundleInteractionReader
          ->_listenedFiberInfo.getAntFiberPointMeshClosestPointIndex();
      _intersectionPointsMeshClosestPointsTexLabels_set.insert(
          _labeledTex0[antFiberPointMeshClosestPoint_index]);
    }
    if (not _intersectionPointsMeshClosestPointsTexLabels_set.empty()) {
      set<short>::const_iterator set_it, set_itb, set_ite;
      set_itb = _intersectionPointsMeshClosestPointsTexLabels_set.begin();
      set_ite = _intersectionPointsMeshClosestPointsTexLabels_set.end();
      
      for (set_it = set_itb; set_it != set_ite; ++ set_it) {
        if ( set_it == set_itb) {
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


    
