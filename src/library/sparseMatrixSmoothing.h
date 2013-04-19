#ifndef CONSTELLATION_SPARSEMATRIXSMOOTHING_H
#define CONSTELLATION_SPARSEMATRIXSMOOTHING_H

#include <constellation/sparseMatrix.h>
#include <constellation/connectivities.h>
#include <aims/mesh/surface.h>

namespace constel
{

  void sparseMatrixDiffusionSmoothing( aims::SparseMatrix & matrix,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const std::vector<size_t> & indices );

  void sparseMatrixDiffusionSmoothing( Connectivities * conn_ptr,
    const AimsTimeSurface<3,Void> & mesh, double connectivityThreshold,
    double distanceThreshold, const std::vector<size_t> & indices );

  void sparseMatrixGaussianSmoothing( aims::SparseMatrix & matrix,
    const AimsSurfaceTriangle & aimsMesh, float distthresh,
    float wthresh = 0.0 );

  void sparseMatrixGaussianSmoothingNormed( aims::SparseMatrix & matrix,
    const AimsSurfaceTriangle & aimsMesh, float distthresh,
    float wthresh = 0.0 );

} // namespace constel

#endif // ifndef CONSTELLATION_SPARSEMATRIXSMOOTHING_H

