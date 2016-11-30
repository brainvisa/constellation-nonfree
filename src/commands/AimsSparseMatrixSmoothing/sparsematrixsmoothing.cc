#include <constellation/sparseMatrixSmoothing.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

int main(int argc, const char** argv) {
  string matFilename;
  Reader<AimsSurfaceTriangle> inMeshAimsR;
  string outmatFilename;
  double sigma;
  double threshold = 0.0001;
  bool gaussian = false;
  Reader<TimeTexture<int32_t> > labeltexR;
  int patch = 0;

  AimsApplication app(
      argc, argv, "Sparse matrix smoothing using heat diffusion or gaussian \n\
      smoothing, with the geometry of a mesh. Can typically be used for \n\
      connectivity matrix smoothing. Smoothing is applied for each line, \n\
      each line is a texture for the mesh.");

  app.addOption(
      matFilename, "-i",
      "input sparse matrix");
  app.addOption(
      inMeshAimsR, "-m",
      "input mesh");
  app.addOption(
      outmatFilename, "-o",
      "output sparse matrix");
  app.addOption(
      sigma, "-s",
      "sigma equivalent of the smoothing");
  app.addOption(
      threshold, "-t",
     "threshold smoothing coefficients under this value (default: 0.0001)",
     true);
  app.addOption(
      labeltexR, "-l",
      "label texture");
  app.addOption(
      patch, "-p",
      "patch index (num in label texture)");
  app.addOption(
      gaussian, "-g",
      "use gaussian smoothing. Default is heat diffusion", true);

  try {
    app.initialize();

    AimsSurfaceTriangle mesh;
    inMeshAimsR.read(mesh);
    rc_ptr<SparseOrDenseMatrix> matrix;
    Reader<SparseOrDenseMatrix> mread(matFilename);
    matrix.reset(mread.read());

    if (gaussian) {
      sparseMatrixGaussianSmoothing(*matrix->asSparse(), mesh, sigma,
                                    threshold);
    } else {
      TimeTexture<int32_t> ptex;
      labeltexR.read(ptex);
      sparseMatrixDiffusionSmoothing(matrix, mesh, threshold, sigma,
                                     ptex, patch);
    }

    Writer<SparseOrDenseMatrix> mwrite(outmatFilename);
    mwrite.write(*matrix);

    return EXIT_SUCCESS;
  }
  catch (user_interruption &) {}

  catch (exception & e) {
    cerr << e.what() << endl;
  }

  return EXIT_FAILURE;
}

