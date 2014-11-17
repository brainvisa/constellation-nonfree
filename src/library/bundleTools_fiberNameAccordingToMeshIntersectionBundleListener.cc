#include <constellation/bundleTools.h>
using namespace aims;
using namespace std;

namespace constel {

  FiberNameAccordingToMeshIntersectionBundleListener::FiberNameAccordingToMeshIntersectionBundleListener(
    BundleInteractionReader &bundleInteractionReader,
    const std::string & fileName )
    : _fileName(fileName)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _file.open(_fileName.c_str(), fstream::out);
    if (_file.is_open())
      cout << "MeshIntersection storing file opened: " << _fileName << endl;
  }

  FiberNameAccordingToMeshIntersectionBundleListener::~FiberNameAccordingToMeshIntersectionBundleListener()
  {
    _file.close();
  }

  void FiberNameAccordingToMeshIntersectionBundleListener::fiberTerminated(
    const BundleProducer &, const BundleInfo &, const FiberInfo & )
  {
    std::set<int> intersectionMeshIdentity_set;
    std::vector<int> fiberMeshIntersectionMeshIdentityVector = _bundleInteractionReader->_listenedFiberInfo.getFiberMeshIntersectionMeshIdentityVector();
    for (uint meshId = 0; meshId < fiberMeshIntersectionMeshIdentityVector.size(); ++meshId)
    {
      intersectionMeshIdentity_set.insert(fiberMeshIntersectionMeshIdentityVector[meshId]);
    }
    if (not intersectionMeshIdentity_set.empty())
    {
      std::set<int>::const_iterator set_it, set_itb, set_ite;
      set_itb = intersectionMeshIdentity_set.begin();
      set_ite = intersectionMeshIdentity_set.end();
      
      for ( set_it = set_itb; set_it != set_ite; ++ set_it)
      {
        if (  set_it == set_itb )
        {
          _file << *set_it << std::flush;
        }
        else
        {
          _file << "_" << *set_it << std::flush;
        }
      }
      _file << std::endl;
    }
    else //intersectionMeshIdentity_set is empty()
    {
      _file << "trash" << std::endl;
    }
  }


} //namespace constel
