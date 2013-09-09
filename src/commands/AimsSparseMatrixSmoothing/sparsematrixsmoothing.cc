#include <constellation/sparseMatrixSmoothing.h>
#include <aims/getopt/getopt2.h>
#include <aims/mesh/texture.h>

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
  Reader<TimeTexture<int32_t> > labeltexR;
  int patch = 0;
  bool iterative_smoothing = false;

  AimsApplication app( argc, argv, "Sparse matrix smoothing using heat diffusion or gaussian smoothing, with the geometry of a mesh. Can typically be used for connectivity matrix smoothing. Smoothing is applied for each line, each line is a texture for the mesh." );
  app.addOption( matFilename, "-i", "input sparse matrix" );
  app.addOption( inMeshAimsR, "-m", "input mesh" );
  app.addOption( outmatFilename, "-o", "output sparse matrix" );
  app.addOption( sigma, "-s", "sigma equivalent of the smoothing" );
  app.addOption( threshold, "-t",
    "threshold smoothing coefficients under this value (default: 0.0001)",
    true );
  app.addOption( labeltexR, "-l", "label texture" );
  app.addOption( patch, "-p", "patch index (num in label texture)" );
  app.addOption( gaussian, "-g",
    "use gaussian smoothing. Default is heat diffusion", true );
  app.addOption( iterative_smoothing, "--iterative_smoothing",
    "use iterative smoothing for the connectivity matrix. By default the diffusion smoothing is implemented using a sparse laplacian coefficient matrix elevated to the power of the number of iterations. It is more efficient, but needs lots of memory for large matrixes. Use this option if your machine lacks memory.", true );

  try
  {
    app.initialize();

    AimsSurfaceTriangle mesh;
    inMeshAimsR.read( mesh );
    rc_ptr<SparseOrDenseMatrix> matrix;
    Reader<SparseOrDenseMatrix> mread( matFilename );
    matrix.reset( mread.read() );

    if( gaussian )
      sparseMatrixGaussianSmoothing( *matrix->asSparse(), mesh, sigma,
                                     threshold );
    else
    {
      TimeTexture<int32_t> ptex;
      labeltexR.read( ptex );
//       vector<size_t> pindices;
//       const vector<int32_t> & labels = ptex[0].data();
//       pindices.reserve( labels.size() ); // arbitrary, should be smaller
//       size_t i, n = labels.size();
//       for( i=0; i!=n; ++i )
//         if( labels[i] == patch )
//           pindices.push_back( i );

      sparseMatrixDiffusionSmoothing( matrix, mesh, threshold, sigma,
                                      ptex, patch, !iterative_smoothing );
    }

    Writer<SparseOrDenseMatrix> mwrite( outmatFilename );
    mwrite.write( *matrix );

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

