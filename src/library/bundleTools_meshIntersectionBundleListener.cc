#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>
#include <constellation/textureAndMeshTools.h>
// includes from CATHIER
#include <cathier/triangle_mesh_geodesic_map.h>

using namespace constel;
using namespace aims;
using namespace std;


namespace constel {
  
  MeshIntersectionBundleListener::MeshIntersectionBundleListener(const AimsSurfaceTriangle & aimsMesh, BundleInteractionReader &bundleInteractionReader, double meshDistanceThreshold, double meshClosestPointMaxDistance, bool verbose ): _aimsMesh(aimsMesh), _meshDistanceThreshold(meshDistanceThreshold), _meshClosestPointMaxDistance(meshClosestPointMaxDistance), _verbose(verbose)
  {
    _bundleInteractionReader = &bundleInteractionReader;
    _meshPolygonsByVertex_Index  = constel::surfacePolygonsIndex( _aimsMesh );
    til::Mesh1 mesh0;
    til::convert(mesh0, _aimsMesh);
    _mesh = addNeighborsToMesh(mesh0);
    // Comuting geomap : neighborhood map
    if (_verbose) std::cout << "meshClosestPoint_minDistance:" << _meshClosestPointMaxDistance <<std::endl;
    if (_verbose) std::cout << "Computing geomap..." << std::flush;
    if (_verbose) std::cout << "step 1..." << std::flush;
    _meshDistanceThresholdNeighborhoodByVertex.reserve(getVertices(_mesh).size());
    for(uint v = 0; v < getVertices(_mesh).size(); ++v)
    {
      _meshDistanceThresholdNeighborhoodByVertex.push_back(QuickMap());
    }
    if (_meshDistanceThreshold > 0)
    {
      til::ghost::GMapStop_AboveThreshold<double> stopGhost(_meshDistanceThreshold);
      boost::shared_ptr<CNeighborhoods> pneighc = til::circular_neighborhoods(getVertices(_mesh), getFaceIndices(_mesh));
      til::Triangle_mesh_geodesic_map<Mesh::VertexCollection, CNeighborhoods, double, til::ghost::GMapStop_AboveThreshold<double>, til::policy::GMap_DefaultStorage_sparse_vect_dbl >
          geomap(getVertices(_mesh), *pneighc, stopGhost);
      std::vector<std::size_t> startPoints(1);
      std::vector<double> dist(1, 0.0);
      std::vector<std::size_t> nneigh(til::size(getVertices(_mesh)));
      {
        for (std::size_t i = 0; i < til::size(getVertices(_mesh)); ++i)
        {
          startPoints[0] = i;
          geomap.init(startPoints, dist);
          geomap.process();
          boost::shared_ptr<til::sparse_vector<double> > tmp = geomap.distanceMap();
          _meshDistanceThresholdNeighborhoodByVertex[i].resize(tmp->getMap().size());
          {
            using namespace til::expr;
            til::detail::loop_xx(castTo(*_1, *_2), _meshDistanceThresholdNeighborhoodByVertex[i], tmp->getMap());
          }
  
          nneigh[i] = _meshDistanceThresholdNeighborhoodByVertex[i].size();
        }
      }
      if (_verbose) std::cout << "OK" << std::endl;
    }
    //KDTREE creation:
    _mesh_kdt_ptr = new KDTree(getVertices(_mesh));
    makeKDTree(getVertices(_mesh), *_mesh_kdt_ptr);
    _mesh_fc_ptr = new til::Find_closest< double, KDTree >(*_mesh_kdt_ptr);
  }
  
  MeshIntersectionBundleListener::~MeshIntersectionBundleListener()
  {
  }

  void MeshIntersectionBundleListener::fiberStarted( const BundleProducer &,
                                  const BundleInfo &,
                                  const FiberInfo & )
  {
    _fiberPointCount = 0;
    _antFiberPoint_ExistingMeshIntersection = false;
    _bundleInteractionReader->_listenedFiberInfo.clearFiberMeshIntersectionInfo();
//     _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(false);
  }

  void MeshIntersectionBundleListener::newFiberPoint( const BundleProducer &,
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
    til::numeric_array<float, 3> fiberPoint_na(fiberPoint[0], fiberPoint[1], fiberPoint[2]);
    //Mesh closest point computing:
    meshClosestPoint_index = (*_mesh_fc_ptr)(fiberPoint_na);
    meshClosestPoint_dist = til::dist2(fiberPoint_na, getVertices(_mesh)[meshClosestPoint_index], til::prec<float>());
    if ( fiberCurvilinearAbscissa > 1 )
    {
      constel::QuickMap * meshPolygonVerticesDistMap_ptr;
      bool intersectMesh;
      if (_meshDistanceThreshold == 0)
      {
        intersectMesh = constel::computeIntersectionPointFiberSegmentAndMesh(_aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint, meshClosestPoint_index, & meshPolygonVerticesDistMap_ptr );
      }
      else //distanceThreshold >0: add neighboorhood to cortexPolygonVerticesDistMap
      {
        intersectMesh = constel::computeIntersectionPointNeighborhoodFiberSegmentAndMesh(_aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint, meshClosestPoint_index, _meshDistanceThresholdNeighborhoodByVertex, & meshPolygonVerticesDistMap_ptr );
      }
      if (intersectMesh)
      {
        //         std::cout << "intersect Cortex" << std::endl;
        _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(*meshPolygonVerticesDistMap_ptr, _meshIdentity);
        delete meshPolygonVerticesDistMap_ptr;
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
          meshPolygonVerticesDistMap_ptr = new QuickMap();
          constel::QuickMap & meshPolygonVerticesDistMap = *meshPolygonVerticesDistMap_ptr;
          //           std::cout <<"Coucou1" << std::endl;
          meshPolygonVerticesDistMap = _meshDistanceThresholdNeighborhoodByVertex[closestPoint_index];
          _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(*meshPolygonVerticesDistMap_ptr, _meshIdentity);
          delete meshPolygonVerticesDistMap_ptr;
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
  
  void MeshIntersectionBundleListener::fiberTerminated( const BundleProducer &, const BundleInfo &, const FiberInfo & )
  {
    if ( _antFiberPoint_ExistingMeshIntersection == false)//if no intersection with cortex, test if dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is the case, add polygones around the meshClosestPoint, according to distanceThreshold, if the previous fiber point has not taken part of an (nearly or not) intersection
    {
//       if( _antFiberPointMeshClosestPoint_dist <= _meshClosestPointMaxDistance and _meshDistanceThreshold!=0 )
      if( _antFiberPointMeshClosestPoint_dist <= _meshClosestPointMaxDistance)
      {
        constel::QuickMap * meshPolygonVerticesDistMap_ptr = new QuickMap();
        constel::QuickMap & meshPolygonVerticesDistMap = *meshPolygonVerticesDistMap_ptr;
        meshPolygonVerticesDistMap = _meshDistanceThresholdNeighborhoodByVertex[_antFiberPointMeshClosestPoint_index];
        _bundleInteractionReader->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(meshPolygonVerticesDistMap, _meshIdentity);
        delete meshPolygonVerticesDistMap_ptr;
        _antFiberPoint_ExistingMeshIntersection = true;
      }
    }
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(_antFiberPoint_ExistingMeshIntersection);
    _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(_antFiberPointMeshClosestPoint_index);
  }

  void MeshIntersectionBundleListener::noMoreBundle( const BundleProducer &)
  {
//     std::cout << "_fiberCount:" << _fiberCount << std::endl;
//     std::cout << "_bundleCortexIntersections nb:" << _bundleCortexIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexNearlyIntersections nb:" << _bundleCortexNearlyIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexConnections nb:" << _bundleCortexConnectionsCount << std::endl;
    delete _mesh_kdt_ptr;
    delete _mesh_fc_ptr;
  }

  void MeshIntersectionBundleListener::setMeshIdentity(int meshIdentity)
  {
    _meshIdentity = meshIdentity;
  }
} // namespace constel
