#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>

using namespace std;
using namespace comist;

namespace constel
{

  MeshIntersectionNoSmoothingFasterBundleListener::MeshIntersectionNoSmoothingFasterBundleListener(const AimsSurfaceTriangle & aimsMesh, AimsData<short> &roisMask, BundleInteractionReader &bundleInteractionReader, double meshClosestPointMaxDistance, bool verbose ): MeshIntersectionNoSmoothingBundleListener(aimsMesh, bundleInteractionReader, meshClosestPointMaxDistance, verbose), _roisMask(roisMask)
  {

  }
  
  MeshIntersectionNoSmoothingFasterBundleListener::~MeshIntersectionNoSmoothingFasterBundleListener()
  {
  }

  void MeshIntersectionNoSmoothingFasterBundleListener::fiberStarted( const comist::BundleProducer & bundleProducer,
                                  const comist::BundleInfo & bundleInfo,
                                  const comist::FiberInfo & fiberPoint)
  {
    MeshIntersectionNoSmoothingBundleListener::fiberStarted(bundleProducer,
                                  bundleInfo,
                                  fiberPoint);
    _antFiberPoint_inRoisMask = false;
  }

  void MeshIntersectionNoSmoothingFasterBundleListener::newFiberPoint( const comist::BundleProducer &,
    const comist::BundleInfo &,
    const comist::FiberInfo &,
    const comist::FiberPoint & fiberPoint)
  {
    ++_fiberPointCount;
    comist::FiberPoint antFiberPoint = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    std::size_t meshClosestPoint_index;
    float meshClosestPoint_dist;
    bool fiberPoint_ExistingMeshIntersection = false;
    float fiberCurvilinearAbscissa = _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa();
    if (_fiberCount <5 and _fiberPointCount <5)
    {
      if (_verbose) std::cout << "fiberCurvilinearAbscissa:" << fiberCurvilinearAbscissa << std::endl;
    }

    //test if the fiberPoint is in the input rois mask (ribbon around the input Mesh, given by _roisMask)
    bool fiberPointInRoisMask = false;
    float x,y,z;
    int indX, indY, indZ;
    short fiberPointInRoisMaskValue;
    x = fiberPoint[0]/_roisMask.sizeX();
    y = fiberPoint[1]/_roisMask.sizeY();
    z = fiberPoint[2]/_roisMask.sizeZ();
    indX = (int) rint(x);
    indY = (int) rint(y);
    indZ = (int) rint(z);
    fiberPointInRoisMaskValue = _roisMask(indX,indY,indZ);
    if ( fiberPointInRoisMaskValue > 0 )
    {
      fiberPointInRoisMask = true;
    }
    //end of the test
    
    if ( fiberPointInRoisMask )
    {
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
        //Compute meshClosestPoint_index and meshClosestPoint_dist of the antFiberPoints if the computation has not been done in the previous step:
        if ( _antFiberPoint_inRoisMask == false)
        {
          til::numeric_array<float, 3> antFiberPoint_na(antFiberPoint[0], antFiberPoint[1], antFiberPoint[2]);
          _antFiberPointMeshClosestPoint_index = (*_mesh_fc_ptr)(antFiberPoint_na);
          _antFiberPointMeshClosestPoint_dist = til::dist2(antFiberPoint_na, getVertices(_mesh)[_antFiberPointMeshClosestPoint_index], til::prec<float>());

        }
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
      
      _antFiberPointMeshClosestPoint_index = meshClosestPoint_index;
      _antFiberPointMeshClosestPoint_dist = meshClosestPoint_dist;
    }
    if (_antFiberPoint_ExistingMeshIntersection)
    {
      _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(_antFiberPoint_ExistingMeshIntersection);
      _bundleInteractionReader->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(_antFiberPointMeshClosestPoint_index);
    }
    _antFiberPoint_inRoisMask = fiberPointInRoisMask;
    _antFiberPoint_ExistingMeshIntersection = fiberPoint_ExistingMeshIntersection;
  }
  
  void MeshIntersectionNoSmoothingFasterBundleListener::fiberTerminated( const comist::BundleProducer &, const comist::BundleInfo &, const comist::FiberInfo & )
  {
    if ( _antFiberPoint_ExistingMeshIntersection == false)//if no intersection with cortex, test if dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is the case, add polygones around the meshClosestPoint, according to distanceThreshold, if the previous fiber point has not taken part of an (nearly or not) intersection
    {
      //Compute meshClosestPoint_index and meshClosestPoint_dist of the antFiberPoints if the computation has not been done in the previous step:
      if ( _antFiberPoint_inRoisMask == false)
      {
        comist::FiberPoint antFiberPoint = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
        til::numeric_array<float, 3> antFiberPoint_na(antFiberPoint[0], antFiberPoint[1], antFiberPoint[2]);
        _antFiberPointMeshClosestPoint_index = (*_mesh_fc_ptr)(antFiberPoint_na);
        _antFiberPointMeshClosestPoint_dist = til::dist2(antFiberPoint_na, getVertices(_mesh)[_antFiberPointMeshClosestPoint_index], til::prec<float>());

      }
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

  void MeshIntersectionNoSmoothingFasterBundleListener::noMoreBundle( const BundleProducer &)
  {
//     std::cout << "_fiberCount:" << _fiberCount << std::endl;
//     std::cout << "_bundleCortexIntersections nb:" << _bundleCortexIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexNearlyIntersections nb:" << _bundleCortexNearlyIntersectionsCount << std::endl;
//     std::cout << "_bundleCortexConnections nb:" << _bundleCortexConnectionsCount << std::endl;
    delete _mesh_kdt_ptr;
    delete _mesh_fc_ptr;
    if (_verbose) std::cout << "fiberCount: " << _fiberCount << std::endl;
  }

  void MeshIntersectionNoSmoothingFasterBundleListener::setMeshIdentity(int meshIdentity)
  {
    _meshIdentity = meshIdentity;
  }
} // namespace constel
