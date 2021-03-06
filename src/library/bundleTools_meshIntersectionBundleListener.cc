#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>
#include <aims/distancemap/meshdistance_d.h>

using namespace constel;
using namespace aims;
using namespace std;


namespace aims
{
  namespace meshdistance
  {

    // implement DistMapMatrixTraits for vector<QuickMap>
    // and instantiate pairwiseDistanceMaps for this type

    template <>
    class DistMapMatrixTraits<QuickMap>
    {
    public:
      static void resize( QuickMap & mat, size_t size1, size_t size2 = 0 )
      {
      }

      static void storeElement( QuickMap & map, size_t i, float value )
      {
        if( value != FLT_MAX )
          map.push_back( make_pair( i, double( value ) ) );
      }
    };

    template
    void
    pairwiseDistanceMaps( const AimsSurface<3,Void> & mesh,
                          std::vector<QuickMap> & distmaps,
                          const Texture<int16_t> & inittex
                            = Texture<int16_t>(),
                          float max_dist = FLT_MAX );
  }
}


namespace constel
{
  
  MeshIntersectionBundleListener::MeshIntersectionBundleListener(
      const AimsSurfaceTriangle &aimsMesh,
      BundleInteractionReader &bundleInteractionReader,
      double meshDistanceThreshold, double meshClosestPointMaxDistance,
      bool verbose) : _meshDistanceThreshold(meshDistanceThreshold),
                      _meshClosestPointMaxDistance(meshClosestPointMaxDistance),
                      _verbose(verbose),
                      _aimsMesh(aimsMesh)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _meshPolygonsByVertex_Index  = constel::surfacePolygonsIndex(_aimsMesh);
    // Comuting geomap : neighborhood map
    if (_verbose)
        cout << "meshClosestPoint_minDistance:" << _meshClosestPointMaxDistance
          << endl;
    if (_verbose) cout << "Computing geomap..." << flush;
    if (_verbose) cout << "step 1..." << flush;

    size_t nv = aimsMesh.vertex().size();

    if( _meshDistanceThreshold > 0 )
      meshdistance::pairwiseDistanceMaps(
        aimsMesh.begin()->second, _meshDistanceThresholdNeighborhoodByVertex,
        _meshDistanceThreshold );
    else
      _meshDistanceThresholdNeighborhoodByVertex.resize( nv );

    //KDTREE creation:
    KDTreeVertices vert = kdt_vertices( _aimsMesh );
    _mesh_kdt_ptr = new KDTree( vert.begin(), vert.end() );
  }
  
  MeshIntersectionBundleListener::~MeshIntersectionBundleListener()
  {}

  void MeshIntersectionBundleListener::fiberStarted(
      const BundleProducer &, const BundleInfo &, const FiberInfo &)
  {
    _fiberPointCount = 0;
    _antFiberPoint_ExistingMeshIntersection = false;
    _bundleInteractionReader
      ->_listenedFiberInfo.clearFiberMeshIntersectionInfo();
  }

  void MeshIntersectionBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint &fiberPoint)
  {
    ++_fiberPointCount;
    FiberPoint antFiberPoint
      = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    size_t meshClosestPoint_index;
    float meshClosestPoint_dist;
    bool fiberPoint_ExistingMeshIntersection = false;
    float fiberCurvilinearAbscissa
      = _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa();
    //Mesh closest point computing:
    meshClosestPoint_index = _mesh_kdt_ptr->find_nearest(
      make_pair( 0U, fiberPoint ) ).first->first;
    meshClosestPoint_dist = dist2(
      fiberPoint, _aimsMesh.vertex()[meshClosestPoint_index] );
    if (fiberCurvilinearAbscissa > 1)
    {
      constel::QuickMap * meshPolygonVerticesDistMap_ptr;
      bool intersectMesh;
      if (_meshDistanceThreshold == 0)
      {
        intersectMesh = constel::computeIntersectionPointFiberSegmentAndMesh(
            _aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint,
            meshClosestPoint_index, & meshPolygonVerticesDistMap_ptr);
      }
      else
      { // >0: add neighboorhood to cortexPolygonVerticesDistMap
        intersectMesh
          = constel::computeIntersectionPointNeighborhoodFiberSegmentAndMesh(
              _aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint,
              fiberPoint, meshClosestPoint_index,
              _meshDistanceThresholdNeighborhoodByVertex,
              &meshPolygonVerticesDistMap_ptr);
      }
      if (intersectMesh)
      {
        _bundleInteractionReader
          ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
              *meshPolygonVerticesDistMap_ptr, _meshIdentity);
        delete meshPolygonVerticesDistMap_ptr;
        fiberPoint_ExistingMeshIntersection =  true;
        _antFiberPoint_ExistingMeshIntersection = true;
      } else if (_antFiberPoint_ExistingMeshIntersection == false
                 and _fiberPointCount == 2)
      {
        // if no intersection with cortex, test if
        // dist(fiberPoint, meshClosestPoint) < dist_min, if it is the case,
        // add polygones around the meshClosestPoint, according to
        // meshDistanceThreshold, if the previous fiber point has not taken
        // part of an (nearly or not) intersection
        float meshClosestPoint_min = min(
            meshClosestPoint_dist, _antFiberPointMeshClosestPoint_dist);
        if (meshClosestPoint_min <= _meshClosestPointMaxDistance)
        {
          size_t closestPoint_index;
          if (meshClosestPoint_dist == meshClosestPoint_min)
          {
            closestPoint_index = meshClosestPoint_index;
            fiberPoint_ExistingMeshIntersection = true;
          }
          else
          { // _antFiberPointMeshClosestPoint_dist==meshClosestPoint_min
            closestPoint_index = _antFiberPointMeshClosestPoint_index;
            _antFiberPoint_ExistingMeshIntersection = true;
          }
          meshPolygonVerticesDistMap_ptr = new QuickMap();
          constel::QuickMap &meshPolygonVerticesDistMap
            = *meshPolygonVerticesDistMap_ptr;
          meshPolygonVerticesDistMap
            = _meshDistanceThresholdNeighborhoodByVertex[closestPoint_index];
          _bundleInteractionReader
            ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
                *meshPolygonVerticesDistMap_ptr, _meshIdentity);
          delete meshPolygonVerticesDistMap_ptr;
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
  
  void MeshIntersectionBundleListener::fiberTerminated(
      const BundleProducer &, const BundleInfo &, const FiberInfo &)
  {
    // if no intersection with cortex, test if 
    // dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is the case,
    // add polygones around the meshClosestPoint, according to
    // distanceThreshold, if the previous fiber point has not taken part of an
    // (nearly or not) intersection
    if (_antFiberPoint_ExistingMeshIntersection == false) {
      if (_antFiberPointMeshClosestPoint_dist
          <= _meshClosestPointMaxDistance) {
        constel::QuickMap * meshPolygonVerticesDistMap_ptr = new QuickMap();
        constel::QuickMap & meshPolygonVerticesDistMap
          = *meshPolygonVerticesDistMap_ptr;
        meshPolygonVerticesDistMap
          = _meshDistanceThresholdNeighborhoodByVertex[
              _antFiberPointMeshClosestPoint_index];
        _bundleInteractionReader
          ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
              meshPolygonVerticesDistMap, _meshIdentity);
        delete meshPolygonVerticesDistMap_ptr;
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

  void MeshIntersectionBundleListener::noMoreBundle(const BundleProducer &) {
    delete _mesh_kdt_ptr;
  }

  void MeshIntersectionBundleListener::setMeshIdentity(int meshIdentity) {
    _meshIdentity = meshIdentity;
  }

} // namespace constel

