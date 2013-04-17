#ifndef CONSTELLATION_CONNMATRIX_H
#define CONSTELLATION_CONNMATRIX_H

#include <constellation/bundleSet.h>
#include <constellation/tildefs.h>
#include <constellation/connectivities.h>

// includes from AIMS
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>


namespace constel
{

  Connectivities * connMatrix( const Fibers & fibers,
                               const AimsSurfaceTriangle & inAimsMesh,
                               float distthresh, float wthresh, Motion motion,
                               bool verbose = false );
  Connectivity * connMatrixSumRows( Connectivities * matrix_ptr,
                                    bool verbose = false );
  Connectivities * connMatrixTargetsRegroup(
    Connectivities * connMatrixToAllMesh,
    const TimeTexture<short> & targetRegionsTexture, int targetRegionsNb,
    bool verbose = false );
  Connectivities * connMatrixRegionExtract(
    Connectivities * allMeshConnMatrix,
    const TimeTexture<short> & seedRegionsTexture, int seedRegionLabel,
    std::size_t seedRegionLabelVertexNb,
    std::vector<std::size_t> ** seedVertexIndex = 0, bool verbose = false );
  Connectivities * connMatrixReducedFromRegion(
    Connectivities * allMeshConnMatrix,
    const TimeTexture<short> & seedRegionsTexture, int seedRegionLabel,
    std::size_t seedRegionLabelVertexNb,
    std::vector<std::size_t> ** seedVertexIndex = 0, bool verbose = false );
//   Connectivities * connMatrixRegionExtract(Connectivities * allMeshConnMatrix, const TimeTexture<short> & seedRegionsTexture, int seedRegionLabel, std::size_t seedRegionLabelVertexNb);
  Connectivities * connMatrixRegionExtractTargetsRegroup(
      Connectivities * allMeshConnMatrix,
      const TimeTexture<short> & seedRegionsTexture, int seedRegionLabel,
      const TimeTexture<short> & targetRegionsTexture, int targetRegionsNb,
      std::size_t seedRegionLabelVertexNb,
      std::vector<std::size_t> ** seedVertexIndex, bool verbose = false );
  void writeAimsFmtConnMatrix( Connectivities * connMatrix_ptr,
                               std::string file_name, bool verbose = false );
  Connectivities * connMatrixNormalize( Connectivities * connMatrix,
                                        bool verbose = false );
  TimeTexture<float> densityTexture(
    Connectivities * allMeshConnMatrix, std::vector<std::size_t> VertexIndex,
    bool verbose = false );
  TimeTexture<float> meshDensityTexture(
    Connectivities * connMatrixToAllMesh_ptr,bool verbose = false );
  TimeTexture<float> * oneTargetDensityTargetsRegroupTexture(
    Connectivity * lineMatrixToTargetRegions_ptr,
    const TimeTexture<short> & targetRegionsTex );
  TimeTexture<float> * connMatrixRow_TO_TimeTexture_FLOAT(
    Connectivity * conn_ptr );
  Connectivity * connMatrixToRois(
    const Fibers & fibers, const AimsSurfaceTriangle & inAimsMesh,
    float distthresh, float wthresh, Motion motion, bool verbose = false );
  Connectivities * connMatrixSeedMeshToTargetMesh(
    const Fibers & fibers, const AimsSurfaceTriangle & aimsMesh1,
    const AimsSurfaceTriangle & aimsMesh2, float distthresh, float wthresh,
    Motion motion, bool verbose = false );
  Connectivities * connMatrixSeedMesh_to_targetMeshTargets_regroup(
    Connectivities * connMatrixSeedMeshToTargetMesh_ptr,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    bool verbose = false );
  Connectivities * connMatrixSeedMeshRegions_to_targetMeshTargets_regroup(
    Connectivities * connMatrixSeedMeshToTargetMesh_ptr,
    const TimeTexture<short> & seedRegionsTexture,
    const TimeTexture<short> & targetRegionsTex, int targetRegionsNb,
    int seedRegionsNb, bool verbose = false );


} // namespace constel

#endif // ifndef CONSTELLATION_CONNMATRIX_H
