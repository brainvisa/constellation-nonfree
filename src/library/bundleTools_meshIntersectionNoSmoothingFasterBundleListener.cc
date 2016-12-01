#include <constellation/bundleTools.h>
#include <constellation/connMatrixTools.h>

using namespace std;
using namespace aims;

namespace constel {

  //---------------------------------------------------------------------------
  MeshIntersectionNoSmoothingFasterBundleListener
    ::MeshIntersectionNoSmoothingFasterBundleListener(
        const AimsSurfaceTriangle &aimsMesh, AimsData<short> &roisMask,
        BundleInteractionReader &bundleInteractionReader,
        double meshClosestPointMaxDistance, bool verbose)
        : MeshIntersectionNoSmoothingBundleListener(
            aimsMesh, bundleInteractionReader, meshClosestPointMaxDistance,
            verbose),
        _roisMask(roisMask) {}

  //---------------------------------------------------------------------------
  MeshIntersectionNoSmoothingFasterBundleListener
    ::~MeshIntersectionNoSmoothingFasterBundleListener() {}

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingFasterBundleListener::fiberStarted(
      const BundleProducer &bundleProducer, const BundleInfo &bundleInfo,
      const FiberInfo &fiberPoint) {
    MeshIntersectionNoSmoothingBundleListener::fiberStarted(bundleProducer,
                                                            bundleInfo,
                                                            fiberPoint);
    _antFiberPoint_inRoisMask = false;
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingFasterBundleListener::newFiberPoint(
      const BundleProducer &, const BundleInfo &, const FiberInfo &,
      const FiberPoint &fiberPoint) {
    ++_fiberPointCount;
    FiberPoint antFiberPoint
      = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
    size_t meshClosestPoint_index;
    float meshClosestPoint_dist;
    bool fiberPoint_ExistingMeshIntersection = false;
    float fiberCurvilinearAbscissa
      = _bundleInteractionReader->_listenedFiberInfo.getCurvilinearAbscissa();
    if (_fiberCount <5 and _fiberPointCount <5) {
      if (_verbose)
        cout << "fiberCurvilinearAbscissa:" << fiberCurvilinearAbscissa
        << endl;
    }

    // test if the fiberPoint is in the input rois mask (ribbon around the
    // input Mesh, given by _roisMask)
    bool fiberPointInRoisMask = false;
    float x, y, z;
    int indX, indY, indZ;
    short fiberPointInRoisMaskValue;
    x = fiberPoint[0]/_roisMask.sizeX();
    y = fiberPoint[1]/_roisMask.sizeY();
    z = fiberPoint[2]/_roisMask.sizeZ();
    indX = (int) rint(x);
    indY = (int) rint(y);
    indZ = (int) rint(z);
    fiberPointInRoisMaskValue = _roisMask(indX,indY,indZ);
    if (fiberPointInRoisMaskValue > 0) {
      fiberPointInRoisMask = true;
    }
    //end of the test
    
    if (fiberPointInRoisMask) {
      til::numeric_array<float, 3> fiberPoint_na(fiberPoint[0],
                                                 fiberPoint[1],
                                                 fiberPoint[2]);
      //Mesh closest point computing:
      meshClosestPoint_index = (*_mesh_fc_ptr)(fiberPoint_na);
      meshClosestPoint_dist = til::dist2(
          fiberPoint_na, getVertices(_mesh)[meshClosestPoint_index],
          til::prec<float>());
      
      if (fiberCurvilinearAbscissa > 1) {
        //Compute meshClosestPoint_index and meshClosestPoint_dist of the
        // antFiberPoints if computation has not been done in the previous step
        if (_antFiberPoint_inRoisMask == false) {
          til::numeric_array<float, 3> antFiberPoint_na(antFiberPoint[0],
                                                        antFiberPoint[1],
                                                        antFiberPoint[2]);
          _antFiberPointMeshClosestPoint_index
            = (*_mesh_fc_ptr)(antFiberPoint_na);
          _antFiberPointMeshClosestPoint_dist = til::dist2(
              antFiberPoint_na,
              getVertices(_mesh)[_antFiberPointMeshClosestPoint_index],
              til::prec<float>());
        }
        constel::QuickMap * meshPolygonVerticesWeightMap_ptr;
        bool intersectMesh;
        intersectMesh = constel::computeIntersectionPointFiberSegmentAndMesh2(
            _aimsMesh, _meshPolygonsByVertex_Index, antFiberPoint, fiberPoint,
            meshClosestPoint_index, & meshPolygonVerticesWeightMap_ptr);
        if (intersectMesh) {
          _bundleInteractionReader
            ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
                *meshPolygonVerticesWeightMap_ptr, _meshIdentity);
          delete meshPolygonVerticesWeightMap_ptr;
          fiberPoint_ExistingMeshIntersection =  true;
          _antFiberPoint_ExistingMeshIntersection = true;
        } else if (_antFiberPoint_ExistingMeshIntersection == false
                    and _fiberPointCount == 2) {
          // if no intersection with cortex,
          // test if dist(fiberPoint, meshClosestPoint) < dist_min, if it is
          // the case, add polygones around the meshClosestPoint, according to
          // meshDistanceThreshold, if the previous fiber point has not taken
          // part of an (nearly or not) intersection
          float meshClosestPoint_min
            = min(meshClosestPoint_dist, _antFiberPointMeshClosestPoint_dist);
          if( meshClosestPoint_min <= _meshClosestPointMaxDistance) {
            size_t closestPoint_index;
            if( meshClosestPoint_dist == meshClosestPoint_min ) {
              closestPoint_index = meshClosestPoint_index;
              fiberPoint_ExistingMeshIntersection = true;
            } else { // _antFiberPointMeshClosestPoint_dist==meshClosestPoint_min
              closestPoint_index = _antFiberPointMeshClosestPoint_index;
              _antFiberPoint_ExistingMeshIntersection = true;
            }
            meshPolygonVerticesWeightMap_ptr = new QuickMap(1);
            constel::QuickMap & meshPolygonVerticesWeightMap
              = *meshPolygonVerticesWeightMap_ptr;
            pair<size_t, double> pair_verticeIndex_WeightClosestToIntersection
              = std::make_pair(closestPoint_index,1.0);
            meshPolygonVerticesWeightMap[0]
              = pair_verticeIndex_WeightClosestToIntersection;
            _bundleInteractionReader
              ->_listenedFiberInfo.pushBackMeshIntersectionNeighbourhood(
                  *meshPolygonVerticesWeightMap_ptr, _meshIdentity);
            delete meshPolygonVerticesWeightMap_ptr;
          }
        }
      }
      _antFiberPointMeshClosestPoint_index = meshClosestPoint_index;
      _antFiberPointMeshClosestPoint_dist = meshClosestPoint_dist;
    }
    if (_antFiberPoint_ExistingMeshIntersection) {
      _bundleInteractionReader
        ->_listenedFiberInfo.setAntFiberPointExistingMeshIntersection(
            _antFiberPoint_ExistingMeshIntersection);
      _bundleInteractionReader
        ->_listenedFiberInfo.setAntFiberPointMeshClosestPointIndex(
            _antFiberPointMeshClosestPoint_index);
    }
    _antFiberPoint_inRoisMask = fiberPointInRoisMask;
    _antFiberPoint_ExistingMeshIntersection = fiberPoint_ExistingMeshIntersection;
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingFasterBundleListener::fiberTerminated(
      const BundleProducer &, const BundleInfo &, const FiberInfo &) {
    // if no intersection with cortex,
    // test if dist(fiberPoint, cortexMeshClosestPoint) < dist_min, if it is
    // the case, add polygones around the meshClosestPoint, according to
    // distanceThreshold, if the previous fiber point has not taken part of an
    // (nearly or not) intersection
    if (_antFiberPoint_ExistingMeshIntersection == false) {
      //Compute meshClosestPoint_index and meshClosestPoint_dist of the
      // antFiberPoints if computation has not been done in the previous step:
      if (_antFiberPoint_inRoisMask == false) {
        FiberPoint antFiberPoint
          = _bundleInteractionReader->_listenedFiberInfo.getAntFiberPoint();
        til::numeric_array<float, 3> antFiberPoint_na(antFiberPoint[0],
                                                      antFiberPoint[1],
                                                      antFiberPoint[2]);
        _antFiberPointMeshClosestPoint_index
          = (*_mesh_fc_ptr)(antFiberPoint_na);
        _antFiberPointMeshClosestPoint_dist
          = til::dist2(antFiberPoint_na,
                       getVertices(_mesh)[_antFiberPointMeshClosestPoint_index],
                       til::prec<float>());

      }
      if (_antFiberPointMeshClosestPoint_dist <= _meshClosestPointMaxDistance) {
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
  void MeshIntersectionNoSmoothingFasterBundleListener::noMoreBundle(
      const BundleProducer &) {
    delete _mesh_kdt_ptr;
    delete _mesh_fc_ptr;
    if (_verbose) std::cout << "fiberCount: " << _fiberCount << std::endl;
  }

  //---------------------------------------------------------------------------
  void MeshIntersectionNoSmoothingFasterBundleListener::setMeshIdentity(
      int meshIdentity) {
    _meshIdentity = meshIdentity;
  }

} // namespace constel
