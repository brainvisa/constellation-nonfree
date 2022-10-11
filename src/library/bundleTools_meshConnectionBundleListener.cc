#include <constellation/bundleTools.h>

using namespace aims;

namespace constel
{

  MeshConnectionBundleListener::MeshConnectionBundleListener(
    BundleInteractionReader &bundleInteractionReader, int meshIdentity,
    bool verbose)
    : _meshIdentity(meshIdentity),
      _verbose(verbose),
      _bundleMeshConnections(new BundleConnections() ),
      _bundleMeshConnectionsLength(new ConnectionsLength() ),
      _bundleMeshConnectionsWeights(std::vector< double >() )
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _bundleMeshConnectionsCount = 0;
  }


  MeshConnectionBundleListener::~MeshConnectionBundleListener() {}


  void MeshConnectionBundleListener::fiberTerminated(
    const BundleProducer &, const BundleInfo &, const FiberInfo &fiberInfo )
  {
    std::vector<constel::QuickMap> fiberIntersectionNeighDistMapVector
      = _bundleInteractionReader
        ->_listenedFiberInfo.getFiberIntersectionNeighDistMapVector();
    std::vector<float> fiberMeshIntersectionCurvilinearAbscissaVector
      = _bundleInteractionReader
        ->_listenedFiberInfo.getFiberMeshIntersectionCurvilinearAbscissaVector();
    std::vector<int> fiberMeshIntersectionMeshIdentityVector
      = _bundleInteractionReader->_listenedFiberInfo.getFiberMeshIntersectionMeshIdentityVector();
    int meshId1, meshId2;
    float l1;
    for (uint intersectionPoint1 = 0;
         intersectionPoint1 < fiberMeshIntersectionCurvilinearAbscissaVector.size();
         ++intersectionPoint1)
    {
      meshId1 = fiberMeshIntersectionMeshIdentityVector[intersectionPoint1];
      l1 = fiberMeshIntersectionCurvilinearAbscissaVector[intersectionPoint1];
      for (uint intersectionPoint2 = intersectionPoint1 + 1;
           intersectionPoint2 < fiberMeshIntersectionCurvilinearAbscissaVector.size();
           ++intersectionPoint2)
      {
        float connection_length
          = fiberMeshIntersectionCurvilinearAbscissaVector[intersectionPoint2]
            - l1;
        if (connection_length != 0)
        {
          meshId2
            = fiberMeshIntersectionMeshIdentityVector[intersectionPoint2];
          Connection connection(2);
          if (meshId1 == meshId2 || meshId1 == _meshIdentity)
          {
            _bundleMeshConnectionsCount++;
            connection[0]
              = fiberIntersectionNeighDistMapVector[intersectionPoint1];
            connection[1]
              = fiberIntersectionNeighDistMapVector[intersectionPoint2];
            _bundleMeshConnections->push_back(connection);
            _bundleMeshConnectionsLength->push_back(connection_length);
            std::cout << "_bundleMeshConnectionsWeights " << fiberInfo.id() << fiberInfo.weight() << std::endl;
            _bundleMeshConnectionsWeights.push_back( fiberInfo.weight() );
          }
        }
      }
    }
  }
} //namespace constel
