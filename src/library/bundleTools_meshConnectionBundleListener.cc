#include <constellation/bundleTools.h>

namespace constel
{

  MeshConnectionBundleListener::MeshConnectionBundleListener(BundleInteractionReader &bundleInteractionReader, int meshIdentity, bool verbose): _meshIdentity(meshIdentity), _bundleMeshConnections(new BundleConnections() ), _bundleMeshConnectionsLength(new ConnectionsLength() ), _verbose(verbose)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _bundleMeshConnectionsCount = 0;
  }
  MeshConnectionBundleListener::~MeshConnectionBundleListener()
  {
  }
  
  void MeshConnectionBundleListener::fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & )
  {
    std::vector<constel::QuickMap > fiberIntersectionNeighDistMapVector = _bundleInteractionReader->_listenedFiberInfo.getFiberIntersectionNeighDistMapVector();
    std::vector<float> fiberMeshIntersectionCurvilinearAbscissaVector = _bundleInteractionReader->_listenedFiberInfo.getFiberMeshIntersectionCurvilinearAbscissaVector();
    std::vector<int> fiberMeshIntersectionMeshIdentityVector =  _bundleInteractionReader->_listenedFiberInfo.getFiberMeshIntersectionMeshIdentityVector();
    int meshId1, meshId2;
    float l1, l2;
    for (uint intersectionPoint1 = 0; intersectionPoint1 < fiberMeshIntersectionCurvilinearAbscissaVector.size(); ++intersectionPoint1)
    {
      meshId1 = fiberMeshIntersectionMeshIdentityVector[intersectionPoint1];
      l1 = fiberMeshIntersectionCurvilinearAbscissaVector[intersectionPoint1];
      for (uint intersectionPoint2 = intersectionPoint1 + 1; intersectionPoint2 < fiberMeshIntersectionCurvilinearAbscissaVector.size(); ++intersectionPoint2)
      {
        float connection_length = fiberMeshIntersectionCurvilinearAbscissaVector[intersectionPoint2] - l1;
        if (connection_length != 0)
        {
//           std::cout << "connection_length:" << connection_length << std::endl;
          meshId2 = fiberMeshIntersectionMeshIdentityVector[intersectionPoint2];
          Connection connection(2);
          if (meshId1 == meshId2 || meshId1 == _meshIdentity)
          {
            _bundleMeshConnectionsCount++;
            connection[0] = fiberIntersectionNeighDistMapVector[intersectionPoint1];
            connection[1] = fiberIntersectionNeighDistMapVector[intersectionPoint2];
            _bundleMeshConnections->push_back(connection);
            _bundleMeshConnectionsLength->push_back(connection_length);
          }
        }
      }
    }
  }
} //namespace constel