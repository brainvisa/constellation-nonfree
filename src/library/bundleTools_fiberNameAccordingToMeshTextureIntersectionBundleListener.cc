#include <constellation/bundleTools.h>
using namespace std;

namespace constel {

  FiberNameAccordingToMeshTextureIntersectionBundleListener::FiberNameAccordingToMeshTextureIntersectionBundleListener(BundleInteractionReader &bundleInteractionReader, const TimeTexture<short> & labeled_tex, std::string fileName ): _fileName(fileName), _labeledTex0(labeled_tex.begin()->second)
  {
    _bundleInteractionReader = &bundleInteractionReader;
//     _labeledTex0 = labeled_tex.begin()->second;
    _file.open(_fileName.c_str(), fstream::out);
    if (_file.is_open()) cout << "MeshIntersection storing file opened: " << _fileName << endl;
  }

  FiberNameAccordingToMeshTextureIntersectionBundleListener::~FiberNameAccordingToMeshTextureIntersectionBundleListener()
  {
    _file.close();
  }

  void FiberNameAccordingToMeshTextureIntersectionBundleListener::fiberStarted( const comist::BundleProducer & bundleProducer,
  const comist::BundleInfo & bundleInfo,
  const comist::FiberInfo & fiberInfo)
  {
    _intersectionPointsMeshClosestPointsTexLabels_set.clear();
  }

  void FiberNameAccordingToMeshTextureIntersectionBundleListener::newFiberPoint( const comist::BundleProducer & bundleProducer, 
  const comist::BundleInfo & bundleInfo,
  const comist::FiberInfo & fiberInfo, 
  const comist::FiberPoint &point )
  {
    if ( (_bundleInteractionReader->_listenedFiberInfo.getAntFiberPointExistingMeshIntersection()) )
    {
      std::size_t antFiberPointMeshClosestPoint_index =  _bundleInteractionReader->_listenedFiberInfo.getAntFiberPointMeshClosestPointIndex();
      _intersectionPointsMeshClosestPointsTexLabels_set.insert(_labeledTex0[antFiberPointMeshClosestPoint_index]);
    }
  }


  void FiberNameAccordingToMeshTextureIntersectionBundleListener::fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & )
  {

    if (_bundleInteractionReader->_listenedFiberInfo.getAntFiberPointExistingMeshIntersection())
    {
      std::size_t antFiberPointMeshClosestPoint_index =  _bundleInteractionReader->_listenedFiberInfo.getAntFiberPointMeshClosestPointIndex();
      _intersectionPointsMeshClosestPointsTexLabels_set.insert(_labeledTex0[antFiberPointMeshClosestPoint_index]);
    }
    if (not _intersectionPointsMeshClosestPointsTexLabels_set.empty())
    {
      std::set<short>::const_iterator set_it, set_itb, set_ite;
      set_itb = _intersectionPointsMeshClosestPointsTexLabels_set.begin();
      set_ite = _intersectionPointsMeshClosestPointsTexLabels_set.end();
      
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


    
