#ifndef CONSTELLATION_CONNMATRIXTOOLS_H
#define CONSTELLATION_CONNMATRIXTOOLS_H

#include <constellation/sparseMatrix.h>
#include <constellation/connectivities.h>
#include <constellation/bundleSet.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

namespace constel
{

  void fillconnMatrix(Connectivities * conn, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, uint connectionLength = 1);

  void fillconnMatrixNoSmoothing(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex);

  void fillconnMatrix(aims::SparseMatrix * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, uint connectionLength= 1);

  void fillconnMatrix(aims::SparseMatrix * conn_ptr, aims::SparseMatrix * conn_ptr2, QuickMap & fiberExtremity1NeighMeshVertex, QuickMap & fiberExtremity2NeighMeshVertex, double connectivityThreshold, double distanceThreshold, std::size_t rowIndex_min, std::size_t rowIndex_max = 0, uint connectionLength= 1);

  void fillNonSymetricConnMatrix(Connectivities * conn_ptr, QuickMap & fiberExtremity1NeighMeshVertex_rows, QuickMap & fiberExtremity2NeighMeshVertex_cols, double connectivityThreshold, double distanceThreshold);

  void fillconnMatrixWithConnections(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold);

//   void fillconnMatrixWithConnections(aims::SparseMatrix * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, std::size_t rowIndex_min = 0, std::size_t rowIndex_max = 0);

  void fillconnMatrixWithConnections(aims::SparseMatrix * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, std::size_t rowIndex_min = 0, std::size_t rowIndex_max = 0, aims::SparseMatrix * conn_ptr2 = 0);

  void fillconnMatrixWithConnectionsPlusLength(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, uint length_min, uint length_max, ConnectionsLength & connectionsLength);
  void fillconnMatrixWithConnectionsPlusLengthWeight(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, uint length_min, uint length_max, ConnectionsLength & connectionsLength);
  void fillconnMatrixWithConnectionsPlusFloatLengthWeight(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold, float length_min, float length_max, ConnectionsFloatLength & connectionsLength);
  void fillNonSymetricConnMatrixWithConnections(Connectivities * conn_ptr, const BundleConnections & connections, double connectivityThreshold, double distanceThreshold);
  template<int D, class T> bool computeIntersectionPointFiberSegmentAndMesh(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** polygonVerticesDistMap = 0 );
  template<int D, class T> bool computeIntersectionPointFiberSegmentAndMesh2(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, QuickMap ** polygonVerticesDistMap = 0 );
  template<int D, class T> bool computeIntersectionPointNeighborhoodFiberSegmentAndMesh(const AimsTimeSurface<D,T> & aimsMesh, const std::vector<std::set<uint> > & polygonsByVertex_Index, Point3df fiberPoint1, Point3df fiberPoint2, uint meshClosestPoint, std::vector<QuickMap> & distanceThresholdNeighborhoodByVertex, QuickMap ** polygonVerticesDistMap_2ptr);

  aims::SparseMatrix* connectivitiesToSparseMatrix(
    const Connectivities & conn );
  Connectivities* sparseMatrixToConnectivities(
    const aims::SparseMatrix & mat );
  void writeConnectivities( const Connectivities & conn,
                            const std::string & filename, bool ascii=false );

  void connMatrixTargetsToTargets(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & targetRegionsTex, std::string filename);

  void connMatrixSeedRegion(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & seedRegionsTex, std::size_t seedRegionLabel, std::string connmatrix_filename, std::string connTexture_filename);

  void connMatrixSeedRegionSmoothed(const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh, Motion motion, const TimeTexture<short> & seedRegionsTex, std::size_t seedRegionLabel, float distthresh, float wthresh, std::string connmatrix_filename, std::string connTexture_filename = "", bool logOption = false);

} // namespace constel

#endif // ifndef CONSTELLATION_CONNMATRIXTOOLS_H
