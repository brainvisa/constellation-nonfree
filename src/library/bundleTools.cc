#include <constellation/bundleTools.h>

using namespace std;
using namespace carto;
using namespace aims;

namespace constel
{

  //-------------//
 //  ListenedFiberInfo  //
//-------------//
  ListenedFiberInfo::ListenedFiberInfo()
  {
    _id = 0;
    _curvilinearAbscissa = 0;
    _antFiberPoint_ExistingMeshIntersection = false;
    _antFiberPointMeshClosestPoint_index = 0;
  }
   ListenedFiberInfo::ListenedFiberInfo( int id ): _id(id)
  {
  }

  void ListenedFiberInfo::setCurvilinearAbscissa(float value)
  {
    _curvilinearAbscissa = value;
  }

 void ListenedFiberInfo::setAntFiberPoint( FiberPoint antFiberPoint)
  {
    _antFiberPoint = antFiberPoint;
  }

  void ListenedFiberInfo::setAntFiberPointExistingMeshIntersection(bool intersect)
  {
    _antFiberPoint_ExistingMeshIntersection = intersect;
  }

  void ListenedFiberInfo::setAntFiberPointMeshClosestPointIndex(std::size_t meshVertex_index)
  {
    _antFiberPointMeshClosestPoint_index = meshVertex_index;
  }
  
  void ListenedFiberInfo::pushBackMeshIntersectionNeighbourhood(constel::QuickMap fiberIntersectionNeighDistMap, int meshIntersectionMeshId)
  {
    _fiberIntersectionNeighDistMapVector.push_back(fiberIntersectionNeighDistMap);
    _fiberMeshIntersectionCurvilinearAbscissaVector.push_back(_curvilinearAbscissa);
    _fiberMeshIntersectionMeshIdentityVector.push_back(meshIntersectionMeshId);
  }
  
  void ListenedFiberInfo::clearFiberMeshIntersectionInfo()
  {
    _fiberIntersectionNeighDistMapVector.clear();
    _fiberMeshIntersectionCurvilinearAbscissaVector.clear();
    _fiberMeshIntersectionMeshIdentityVector.clear();
  }

//   void ListenedFiberInfo::savingFiberMeshIntersectionPoints(std::fstream file)
//   {
//     for (unsigned intersectionPoint = 0; intersectionPoint < _fiberMeshIntersectionCurvilinearAbscissaVector.size(); ++intersectionPoint)
//     {
//       for (unsigned polygon_neigh = 0; polygon_neigh< _fiberIntersectionNeighDistMapVector[intersectionPoint].size(); ++ polygon_neigh)
//       {
//         file << _fiberMeshIntersectionMeshNameVector[intersectionPoint] << " : " << _fiberIntersectionNeighDistMapVector[intersectionPoint][polygon_neigh].first << endl;
//       }
//     }
//   }

  //-------------//
 //  BundleInteractionReader  //
//-------------//
  
// 
  BundleInteractionReader::BundleInteractionReader( const std::string &fileName ): BundleReader( fileName ), _listenedFiberInfo( ListenedFiberInfo() ), _meshIntersectionBundleListener_nb(0)
  {
//     BundleReader( fileName );
    std::cout << "BundleInteractionReader::BundleInteractionReader" << std::endl;
    _meshIntersectionBundleListener_nb = 0;
  }
  BundleInteractionReader::~BundleInteractionReader()
  {
  }
//   void BundleInteractionReader::addBundleListener( MeshIntersectionBundleListener & meshIntersectionBundleListener)
//   {
//     std::cout << "BundleInteractionReader::addBundleListener" << std::endl;
//     std::cout << "_meshIntersectionBundleListener_nb:" << _meshIntersectionBundleListener_nb << std::endl;
//     meshIntersectionBundleListener.setMeshIdentity(_meshIntersectionBundleListener_nb);
//     ++_meshIntersectionBundleListener_nb;
//     BundleProducer::addBundleListener(meshIntersectionBundleListener);

//   }


  //-------------//
 //  MemAntBundleListener  //
//-------------//
  MemAntBundleListener::MemAntBundleListener(BundleInteractionReader &bundleInteractionReader )
  {
    _bundleInteractionReader = &bundleInteractionReader;
  }
  
  MemAntBundleListener::~MemAntBundleListener()
  {
  }

  void MemAntBundleListener::newFiberPoint( const BundleProducer &,
    const BundleInfo &,
    const FiberInfo &,
    const FiberPoint & fiberPoint)
  {
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPoint(fiberPoint);
//     std::cout << "MemAntBundleListener::newFiberPoint: antFiberPoint:" << fiberPoint[0] << ", " << fiberPoint[1] << ", " << fiberPoint[2] << std::endl;
  }

  //-------------//
 //  AfficheAntFiberPointBundleListener  //
//-------------//

  AfficheAntFiberPointBundleListener::AfficheAntFiberPointBundleListener(BundleInteractionReader &bundleInteractionReader )
  {
    _bundleInteractionReader = &bundleInteractionReader;
  }
  void AfficheAntFiberPointBundleListener::newFiberPoint( const BundleProducer &,
    const BundleInfo &,
    const FiberInfo &,
    const FiberPoint & fiberPoint)
  {
    FiberPoint coucouPoint = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    std::cout << "AfficheAntFiberPointBundleListener::newFiberPoint: coucouPoint:" << coucouPoint[0] << ", " << coucouPoint[1] << ", " << coucouPoint[2] << std::endl;
  }
  AfficheAntFiberPointBundleListener::~AfficheAntFiberPointBundleListener()
  {
  }

  //-------------//
 //  CurvilinearAbscissaBundleListener  //
//-------------//

  CurvilinearAbscissaBundleListener::CurvilinearAbscissaBundleListener(BundleInteractionReader &bundleInteractionReader )
  {
    _bundleInteractionReader = &bundleInteractionReader;
  }
  
  void CurvilinearAbscissaBundleListener::fiberStarted( const BundleProducer &, const BundleInfo &, const FiberInfo & )
  {
    _fiberPointCount = 0;
    _fiberLength = 0.0;
    _bundleInteractionReader->_listenedFiberInfo.setCurvilinearAbscissa(0.0);
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
//     std::cout << "CurvilinearAbscissaBundleListener::fiberStarted: curv abscissa:" << _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa() << std::endl;
  }

  void CurvilinearAbscissaBundleListener::newFiberPoint( const BundleProducer &,
    const BundleInfo &,
    const FiberInfo &,
    const FiberPoint & fiberPoint)
  {
    ++_fiberPointCount;
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
    
    if ( _fiberPointCount > 1 )
    {
      til::numeric_array<float, 3> fiberPoint_na(fiberPoint[0], fiberPoint[1], fiberPoint[2]);
      FiberPoint antFiberPoint = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
      til::numeric_array<float, 3> antFiberPoint_na(antFiberPoint[0], antFiberPoint[1], antFiberPoint[2]);
      _fiberLength += sqrt(til::dist2(fiberPoint_na, antFiberPoint_na, til::prec<float>()));
      _bundleInteractionReader->_listenedFiberInfo.setCurvilinearAbscissa(_fiberLength);
//       if ( _fiberPointCount < 10 )
//       {
//       std::cout << "newFiberPoint curv abscissa:" << _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa() << std::endl;
//       }
    }
  }
  CurvilinearAbscissaBundleListener::~CurvilinearAbscissaBundleListener()
  {
  }
}//namespace constel

