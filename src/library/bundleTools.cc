#include <constellation/bundleTools.h>

using namespace std;
using namespace carto;
using namespace aims;

namespace constel {

  //---------------------
  //  ListenedFiberInfo  
  //---------------------
  ListenedFiberInfo::ListenedFiberInfo() {
    _id = 0;
    _curvilinearAbscissa = 0;
    _antFiberPoint_ExistingMeshIntersection = false;
    _antFiberPointMeshClosestPoint_index = 0;
  }

  ListenedFiberInfo::ListenedFiberInfo(int id): _id(id) {}

  void ListenedFiberInfo::setCurvilinearAbscissa(float value) {
    _curvilinearAbscissa = value;
  }

  void ListenedFiberInfo::setAntFiberPoint(FiberPoint antFiberPoint) {
    _antFiberPoint = antFiberPoint;
  }

  void ListenedFiberInfo::setAntFiberPointExistingMeshIntersection(
      bool intersect) {
    _antFiberPoint_ExistingMeshIntersection = intersect;
  }

  void ListenedFiberInfo::setAntFiberPointMeshClosestPointIndex(
      std::size_t meshVertex_index) {
    _antFiberPointMeshClosestPoint_index = meshVertex_index;
  }
  
  void ListenedFiberInfo::pushBackMeshIntersectionNeighbourhood(
      const constel::QuickMap & fiberIntersectionNeighDistMap,
      int meshIntersectionMeshId)
  {
    _fiberIntersectionNeighDistMapVector.push_back(
        fiberIntersectionNeighDistMap);
    _fiberMeshIntersectionCurvilinearAbscissaVector.push_back(
        _curvilinearAbscissa);
    _fiberMeshIntersectionMeshIdentityVector.push_back(
        meshIntersectionMeshId);
  }
  
  void ListenedFiberInfo::clearFiberMeshIntersectionInfo() {
    _fiberIntersectionNeighDistMapVector.clear();
    _fiberMeshIntersectionCurvilinearAbscissaVector.clear();
    _fiberMeshIntersectionMeshIdentityVector.clear();
  }


  //---------------------------
  //  BundleInteractionReader  
  //---------------------------
  BundleInteractionReader::BundleInteractionReader(const std::string &fileName)
      : BundleReader(fileName), _listenedFiberInfo(ListenedFiberInfo())
  {
  }

  BundleInteractionReader::~BundleInteractionReader() {}


  //------------------------
  //  MemAntBundleListener  
  //------------------------
  MemAntBundleListener::MemAntBundleListener(
      BundleInteractionReader &bundleInteractionReader)
  {
    _bundleInteractionReader = &bundleInteractionReader;
  }
  
  MemAntBundleListener::~MemAntBundleListener() {}

  void MemAntBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint & fiberPoint)
  {
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPoint(fiberPoint);
  }


  //--------------------------------------
  //  AfficheAntFiberPointBundleListener  
  //--------------------------------------
  AfficheAntFiberPointBundleListener::AfficheAntFiberPointBundleListener(
      BundleInteractionReader &bundleInteractionReader) {
    _bundleInteractionReader = &bundleInteractionReader;
  }

  void AfficheAntFiberPointBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint & /* fiberPoint */) {
    FiberPoint coucouPoint =
        _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
  }

  AfficheAntFiberPointBundleListener::~AfficheAntFiberPointBundleListener() {}


  //-------------------------------------
  //  CurvilinearAbscissaBundleListener  
  //-------------------------------------
  CurvilinearAbscissaBundleListener::CurvilinearAbscissaBundleListener(
      BundleInteractionReader &bundleInteractionReader)
  {
    _bundleInteractionReader = &bundleInteractionReader;
  }

  void CurvilinearAbscissaBundleListener::fiberStarted(
      const BundleProducer &, const BundleInfo &, const FiberInfo &)
  {
    _fiberPointCount = 0;
    _fiberLength = 0.0;
    _bundleInteractionReader->_listenedFiberInfo.setCurvilinearAbscissa(0.0);
    _bundleInteractionReader->
        _listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
  }

  void CurvilinearAbscissaBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint & fiberPoint)
  {
    ++_fiberPointCount;
    _bundleInteractionReader->
        _listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
    
    if (_fiberPointCount > 1)
    {
      FiberPoint antFiberPoint =
          _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
      _fiberLength += sqrt(dist2( fiberPoint, antFiberPoint ));
      _bundleInteractionReader->_listenedFiberInfo.setCurvilinearAbscissa(
          _fiberLength);
    }
  }

  CurvilinearAbscissaBundleListener::~CurvilinearAbscissaBundleListener() {}
}//namespace constel

