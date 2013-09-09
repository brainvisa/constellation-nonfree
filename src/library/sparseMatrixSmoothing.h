#ifndef CONSTELLATION_SPARSEMATRIXSMOOTHING_H
#define CONSTELLATION_SPARSEMATRIXSMOOTHING_H

#include <aims/sparsematrix/sparseMatrix.h>
#include <constellation/connectivities.h>
#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

namespace constel
{

  void sparseMatrixDiffusionSmoothing( aims::SparseMatrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const TimeTexture<int32_t> & patches,
    int32_t patch, bool use_matpow = true );

  void sparseMatrixDiffusionSmoothing( aims::SparseMatrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const TimeTexture<int16_t> & patches,
    int32_t patch, bool use_matpow = true );

  void sparseMatrixDiffusionSmoothing( Connectivities * conn_ptr,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const TimeTexture<int32_t> & patches,
    int32_t patch, bool use_matpow = true );

  void sparseMatrixDiffusionSmoothing( Connectivities * conn_ptr,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const TimeTexture<int16_t> & patches,
    int32_t patch, bool use_matpow = true );

  void sparseMatrixGaussianSmoothing( aims::SparseMatrix & matrix,
    const AimsSurfaceTriangle & aimsMesh, float distthresh,
    float wthresh = 0.0 );

  void sparseMatrixGaussianSmoothingNormed( aims::SparseMatrix & matrix,
    const AimsSurfaceTriangle & aimsMesh, float distthresh,
    float wthresh = 0.0 );

} // namespace constel

#endif // ifndef CONSTELLATION_SPARSEMATRIXSMOOTHING_H

