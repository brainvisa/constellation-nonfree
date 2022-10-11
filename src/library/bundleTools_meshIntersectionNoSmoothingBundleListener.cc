#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>

using namespace constel;
using namespace aims;
using namespace std;

namespace constel {

  //---------------------------------------------------------------------------
  MeshIntersectionNoSmoothingBundleListener
    ::MeshIntersectionNoSmoothingBundleListener(
        const AimsSurfaceTriangle & aimsMesh,
        BundleInteractionReader &bundleInteractionReader,
        double meshClosestPointMaxDistance, bool verbose)
        : _meshClosestPointMaxDistance(meshClosestPointMaxDistance),
          _verbose(verbose),
          _aimsMesh(aimsMesh)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _meshPolygonsByVertex_Index  = constel::surfacePolygonsIndex(_aimsMesh);
    //KDTREE creation:
    KDTreeVertices vert = kdt_vertices( _aimsMesh );
    _mesh_kdt_ptr = new KDTree( vert.begin(), vert.end() );

    _fiberCount = 0;
  }

  //---------------------------------------------------------------------------
  MeshIntersectionNoSmoothingBundleListener
    ::~MeshIntersectionNoSmoothingBundleListener() {}

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingBundleListener::fiberStarted(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    _fiberPointCount = 0;
    _fiberCount += 1;
    _antFiberPoint_ExistingMeshIntersection = false;
    _bundleInteractionReader
      ->_listenedFiberInfo.clearFiberMeshIntersectionInfo();
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint & fiberPoint)
  {
    ++_fiberPointCount;
    FiberPoint antFiberPoint
      = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    std::size_t meshClosestPoint_index;
    float meshClosestPoint_dist;
    bool fiberPoint_ExistingMeshIntersection = false;
    float fiberCurvilinearAbscissa
      = _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa();
    if (_fiberCount < 5 and _fiberPointCount < 5) {
      if (_verbose)
        cout << "fiberCurvilinearAbscissa:" << fiberCurvilinearAbscissa << endl;
    }
    //Mesh closest point computing:
    meshClosestPoint_index = _mesh_kdt_ptr->find_nearest(
      make_pair( 0U, fiberPoint ) ).first->first;
    meshClosestPoint_dist = dist2(
      fiberPoint, _aimsMesh.vertex()[meshClosestPoint_index] );
    
    if (fiberCurvilinearAbscissa > 1)
    {
      constel::QuickMap * meshPolygonVerticesWeightMap_ptr;
      bool intersectMesh;
      intersectMesh = constel::computeIntersectionPointFiberSegmentAndMesh2(
          _aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint,
          meshClosestPoint_index, & meshPolygonVerticesWeightMap_ptr);
      if (intersectMesh)
      {
        _bundleInteractionReader
          ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
              *meshPolygonVerticesWeightMap_ptr, _meshIdentity);
        delete meshPolygonVerticesWeightMap_ptr;
        fiberPoint_ExistingMeshIntersection =  true;
        _antFiberPoint_ExistingMeshIntersection = true;
      }
      else if (_antFiberPoint_ExistingMeshIntersection == false
                && _fiberPointCount == 2)
      {
        // if no intersection with cortex,
        // test if dist(fiberPoint, meshClosestPoint) < dist_min, if it is the
        // case, add polygones around the meshClosestPoint, according to
        // meshDistanceThreshold, if the previous fiber point has not taken
        // part of an (nearly or not) intersection
        float meshClosestPoint_min
          = min(meshClosestPoint_dist, _antFiberPointMeshClosestPoint_dist);
        if (meshClosestPoint_min <= _meshClosestPointMaxDistance)
        {
          size_t closestPoint_index;
          if ( meshClosestPoint_dist == meshClosestPoint_min)
          {
            closestPoint_index = meshClosestPoint_index;
            fiberPoint_ExistingMeshIntersection = true;
          }
          else
          { //_antFiberPointMeshClosestPoint_dist==meshClosestPoint_min
            closestPoint_index = _antFiberPointMeshClosestPoint_index;
            _antFiberPoint_ExistingMeshIntersection = true;
          }
          meshPolygonVerticesWeightMap_ptr = new QuickMap(1);
          constel::QuickMap & meshPolygonVerticesWeightMap
            = *meshPolygonVerticesWeightMap_ptr;
          pair<size_t, double> pair_verticeIndex_WeightClosestToIntersection
            = make_pair(closestPoint_index, 1.0);
          meshPolygonVerticesWeightMap[0]
            = pair_verticeIndex_WeightClosestToIntersection;
          _bundleInteractionReader
            ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
                *meshPolygonVerticesWeightMap_ptr, _meshIdentity);
          delete meshPolygonVerticesWeightMap_ptr;
        }
      }
    }
    if (_antFiberPoint_ExistingMeshIntersection)
    {
      _bundleInteractionReader
        ->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(
            _antFiberPoint_ExistingMeshIntersection);
      _bundleInteractionReader
        ->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(
            _antFiberPointMeshClosestPoint_index);
    }
    _antFiberPointMeshClosestPoint_index = meshClosestPoint_index;
    _antFiberPointMeshClosestPoint_dist = meshClosestPoint_dist;
    _antFiberPoint_ExistingMeshIntersection
      = fiberPoint_ExistingMeshIntersection;
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingBundleListener::fiberTerminated(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    if (_antFiberPoint_ExistingMeshIntersection == false)//if no intersection with cortex, test if dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is the case, add polygones around the meshClosestPoint, according to distanceThreshold, if the previous fiber point has not taken part of an (nearly or not) intersection
    {
      if (_antFiberPointMeshClosestPoint_dist
          <= _meshClosestPointMaxDistance)
      {
        constel::QuickMap * meshPolygonVerticesWeightMap_ptr = new QuickMap(1);
        constel::QuickMap & meshPolygonVerticesWeightMap
          = *meshPolygonVerticesWeightMap_ptr;
        pair<size_t, double> pair_verticeIndex_WeightClosestToIntersection
          = make_pair(_antFiberPointMeshClosestPoint_index, 1.0);
        meshPolygonVerticesWeightMap[0]
          = pair_verticeIndex_WeightClosestToIntersection;
        _bundleInteractionReader
          ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
              meshPolygonVerticesWeightMap, _meshIdentity);
        delete meshPolygonVerticesWeightMap_ptr;
        _antFiberPoint_ExistingMeshIntersection = true;
      }
    }
    _bundleInteractionReader
      ->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(
          _antFiberPoint_ExistingMeshIntersection);
    _bundleInteractionReader
      ->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(
          _antFiberPointMeshClosestPoint_index);
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingBundleListener::noMoreBundle(
      const BundleProducer &) {
    delete _mesh_kdt_ptr;
    if (_verbose) cout << "fiberCount: " << _fiberCount << endl;
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingBundleListener::setMeshIdentity(
      int meshIdentity) {
    _meshIdentity = meshIdentity;
  }

} // namespace constel
