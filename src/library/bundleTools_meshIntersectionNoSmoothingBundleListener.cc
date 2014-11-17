#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>

using namespace constel;
using namespace aims;
using namespace std;

namespace constel
{

  MeshIntersectionNoSmoothingBundleListener::MeshIntersectionNoSmoothingBundleListener(
    const AimsSurfaceTriangle & aimsMesh,
    BundleInteractionReader &bundleInteractionReader,
    double meshClosestPointMaxDistance, bool verbose )
    : _aimsMesh(aimsMesh),
    _meshClosestPointMaxDistance(meshClosestPointMaxDistance),
    _verbose(verbose)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _meshPolygonsByVertex_Index  = constel::surfacePolygonsIndex( _aimsMesh );
    til::Mesh1 mesh0;
    til::convert(mesh0, _aimsMesh);
    _mesh = addNeighborsToMesh(mesh0);
    //KDTREE creation:
    _mesh_kdt_ptr = new KDTree(getVertices(_mesh));
    makeKDTree(getVertices(_mesh), *_mesh_kdt_ptr);
    _mesh_fc_ptr = new til::Find_closest< double, KDTree >(*_mesh_kdt_ptr);
    _fiberCount = 0;
  }
  
  MeshIntersectionNoSmoothingBundleListener::~MeshIntersectionNoSmoothingBundleListener()
  {
  }

  void MeshIntersectionNoSmoothingBundleListener::fiberStarted(
    const BundleProducer &,
    const BundleInfo &,
    const FiberInfo & )
  {
    _fiberPointCount = 0;
    _fiberCount += 1;
    _antFiberPoint_ExistingMeshIntersection = false;
    _bundleInteractionReader->_listenedFiberInfo.clearFiberMeshIntersectionInfo();
//     _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
  }

  void MeshIntersectionNoSmoothingBundleListener::newFiberPoint(
    const BundleProducer &,
    const BundleInfo &,
    const FiberInfo &,
    const FiberPoint & fiberPoint)
  {
    ++_fiberPointCount;
    FiberPoint antFiberPoint = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    std::size_t meshClosestPoint_index;
    float meshClosestPoint_dist;
    bool fiberPoint_ExistingMeshIntersection = false;
    float fiberCurvilinearAbscissa = _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa();
    if (_fiberCount <5 and _fiberPointCount <5)
    {
      if (_verbose) std::cout << "fiberCurvilinearAbscissa:" << fiberCurvilinearAbscissa << std::endl;
    }
    til::numeric_array<float, 3> fiberPoint_na(fiberPoint[0], fiberPoint[1], fiberPoint[2]);
    //Mesh closest point computing:
//     std::cout << "Look for mesh closest point index..." << std::flush;
    meshClosestPoint_index = (*_mesh_fc_ptr)(fiberPoint_na);
//     std::cout << "Done." << std::endl;
//     std::cout << "Compute distance between fiber point and mesh closest point ..." << std::flush;
    meshClosestPoint_dist = til::dist2(fiberPoint_na, getVertices(_mesh)[meshClosestPoint_index], til::prec<float>());
//     std::cout << "Done. Distance : " << meshClosestPoint_dist << std::endl;
    
    if ( fiberCurvilinearAbscissa > 1 )
    {
      constel::QuickMap * meshPolygonVerticesWeightMap_ptr;
      bool intersectMesh;
//       std::cout << "Compute intersection point between fiber point and polygons around the mesh closest point ..." << std::flush;
      intersectMesh = constel::computeIntersectionPointFiberSegmentAndMesh2(_aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint, meshClosestPoint_index, & meshPolygonVerticesWeightMap_ptr );
//       std::cout << "Done." << std::endl;
      if (intersectMesh)
      {
        //         std::cout << "intersect Cortex" << std::endl;
        _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(*meshPolygonVerticesWeightMap_ptr, _meshIdentity);
        delete meshPolygonVerticesWeightMap_ptr;
        fiberPoint_ExistingMeshIntersection =  true;
        _antFiberPoint_ExistingMeshIntersection = true;
      }
      else if ( _antFiberPoint_ExistingMeshIntersection == false and _fiberPointCount == 2)//if no intersection with cortex, test if dist(fiberPoint, meshClosestPoint) < dist_min, if it is the case, add polygones around the meshClosestPoint, according to meshDistanceThreshold, if the previous fiber point has not taken part of an (nearly or not) intersection
      {
        //         std::cout << "cortexMeshClosestPoint_dist:" << cortexMeshClosestPoint_dist << ", ant_cortexMeshClosestPoint_dist" << _antCortexMeshClosestPoint_dist << std::endl;
        float meshClosestPoint_min = min(meshClosestPoint_dist, _antFiberPointMeshClosestPoint_dist);
        //         std::cout << "min cortexMeshClosestPoint_dist, ant_cortexMeshClosestPoint_dist:" << meshClosestPoint_min << std::endl;
//         if( meshClosestPoint_min <= _meshClosestPointMaxDistance and _meshDistanceThreshold!=0 )
        if( meshClosestPoint_min <= _meshClosestPointMaxDistance)
        {
          size_t closestPoint_index;
          if( meshClosestPoint_dist == meshClosestPoint_min )
          {
            closestPoint_index = meshClosestPoint_index;
            fiberPoint_ExistingMeshIntersection = true;
          }
          else // _antFiberPointMeshClosestPoint_dist == meshClosestPoint_min
          {
            closestPoint_index = _antFiberPointMeshClosestPoint_index;
            _antFiberPoint_ExistingMeshIntersection = true;
          }
          meshPolygonVerticesWeightMap_ptr = new QuickMap(1);
          constel::QuickMap & meshPolygonVerticesWeightMap = *meshPolygonVerticesWeightMap_ptr;
          //           std::cout <<"Coucou1" << std::endl;
          std::pair<std::size_t, double > pair_verticeIndex_WeightClosestToIntersection = std::make_pair(closestPoint_index,1.0);
          meshPolygonVerticesWeightMap[0] = pair_verticeIndex_WeightClosestToIntersection;
          _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(*meshPolygonVerticesWeightMap_ptr, _meshIdentity);
          delete meshPolygonVerticesWeightMap_ptr;
        }
      }
    }
    if (_antFiberPoint_ExistingMeshIntersection)
    {
      _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(_antFiberPoint_ExistingMeshIntersection);
      _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(_antFiberPointMeshClosestPoint_index);
    }
    _antFiberPointMeshClosestPoint_index = meshClosestPoint_index;
    _antFiberPointMeshClosestPoint_dist = meshClosestPoint_dist;
    _antFiberPoint_ExistingMeshIntersection = fiberPoint_ExistingMeshIntersection;
  }
  
  void MeshIntersectionNoSmoothingBundleListener::fiberTerminated(
    const BundleProducer &, const BundleInfo &, const FiberInfo & )
  {
    if ( _antFiberPoint_ExistingMeshIntersection == false)//if no intersection with cortex, test if dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is the case, add polygones around the meshClosestPoint, according to distanceThreshold, if the previous fiber point has not taken part of an (nearly or not) intersection
    {
//       if( _antFiberPointMeshClosestPoint_dist <= _meshClosestPointMaxDistance and _meshDistanceThreshold!=0 )
      if( _antFiberPointMeshClosestPoint_dist <= _meshClosestPointMaxDistance)
      {
        constel::QuickMap * meshPolygonVerticesWeightMap_ptr = new QuickMap(1);
        constel::QuickMap & meshPolygonVerticesWeightMap = *meshPolygonVerticesWeightMap_ptr;
        std::pair<std::size_t, double > pair_verticeIndex_WeightClosestToIntersection = std::make_pair(_antFiberPointMeshClosestPoint_index,1.0);
        meshPolygonVerticesWeightMap[0] = pair_verticeIndex_WeightClosestToIntersection;
        _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(meshPolygonVerticesWeightMap, _meshIdentity);
        delete meshPolygonVerticesWeightMap_ptr;
        _antFiberPoint_ExistingMeshIntersection = true;
      }
    }
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(_antFiberPoint_ExistingMeshIntersection);
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(_antFiberPointMeshClosestPoint_index);
  }

  void MeshIntersectionNoSmoothingBundleListener::noMoreBundle( const BundleProducer &)
  {
//     std::cout << "_fiberCount:" << _fiberCount << std::endl;
//     std::cout << "_bundleCortexIntersections nb:" << _bundleCortexIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexNearlyIntersections nb:" << _bundleCortexNearlyIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexConnections nb:" << _bundleCortexConnectionsCount << std::endl;
    delete _mesh_kdt_ptr;
    delete _mesh_fc_ptr;
    if (_verbose) std::cout << "fiberCount: " << _fiberCount << std::endl;
  }

  void MeshIntersectionNoSmoothingBundleListener::setMeshIdentity(int meshIdentity)
  {
    _meshIdentity = meshIdentity;
  }
} // namespace constel
