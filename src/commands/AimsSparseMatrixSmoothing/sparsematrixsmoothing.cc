#include <constellation/sparseMatrixSmoothing.h>
#include <aims/getopt/getopt2.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

int main( int argc, const char** argv )
{
  string matFilename;
  Reader<AimsSurfaceTriangle> inMeshAimsR;
  string outmatFilename;
  double sigma;
  double threshold = 0.0001;
  bool gaussian = false;

  AimsApplication app( argc, argv, "Sparse matrix smoothing using heat diffusion or gaussian smoothing, with the geometry of a mesh. Can typically be used for connectivity matrix smoothing. Smoothing is applied for each line, each line is a texture for the mesh." );
  app.addOption( matFilename, "-i", "input sparse matrix" );
  app.addOption( inMeshAimsR, "-m", "input mesh" );
  app.addOption( outmatFilename, "-o", "output sparse matrix" );
  app.addOption( sigma, "-s", "sigma equivalent of the smoothing" );
  app.addOption( threshold, "-t",
    "threshold smoothing coefficients under this value (default: 0.0001)",
    true );
  app.addOption( gaussian, "-g",
    "use gaussian smoothing. Default is heat diffusion", true );

  try
  {
    app.initialize();

    AimsSurfaceTriangle mesh;
    inMeshAimsR.read( mesh );
    SparseMatrix matrix;
    matrix.read( matFilename );

    if( gaussian )
      sparseMatrixGaussianSmoothing( matrix, mesh, sigma, threshold );
    else
      sparseMatrixDiffusionSmoothing( matrix, mesh, threshold, sigma );
    matrix.write( outmatFilename );

    return EXIT_SUCCESS;
  }
  catch( user_interruption & )
  {
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
  return EXIT_FAILURE;
}

