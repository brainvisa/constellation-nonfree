#include <constellation/connMatrix.h>
#include <aims/sparsematrix/sparseMatrix.h>
#include <aims/mesh/texturetools.h>
#include <aims/getopt/getopt2.h>
#include <aims/sparsematrix/sparseordensematrix.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace comist;
using namespace constel;


int main( int argc, const char** argv )
{
  try
  {
    typedef AimsData<float> Matrix;

    Reader< AimsSurfaceTriangle > inMeshAimsR;
    Reader< TimeTexture<short> > seedRegionsTexR;
    Reader< TimeTexture<short> > targetRegionsTexR;
    AimsSurfaceTriangle inAimsMesh;
    TimeTexture<short> seedRegionsTex;
    TimeTexture<short> targetRegionsTex;
    size_t targetRegionsNb = 0;
    size_t seedRegionLabel = 0;
    size_t seedRegionsNb = 0;
    string connTextureFileName;
    string connMatrixFileName = "";
    bool verbose = false;
    bool normalize = false;
    string connectivityTextureType = "oneSeedRegion_to_targets";
    string connLineTargetsDensityFileName = "";
    string inConnMatrixFileName = "";

    AimsApplication app( argc, argv, "computes and writes connection density texture(s).\n\
    Modes supported: (1) oneSeedRegion_to_targets (2) seedVertex_to_targets" );
    app.addOption( inMeshAimsR, "-mesh", "input mesh" );
    app.addOption( inConnMatrixFileName, "-connmatrixfile",
                   "input conn matrix filename." );
    app.addOption( seedRegionsTexR, "-seedregionstex",
                   "input region texture : study a region in particular." );
    app.addOption( connTextureFileName, "-outconntex",
                   "output mean connectivity texture file name.", true );
    app.addOption( seedRegionLabel, "-seedlabel",
                   "input seed region label.\n\
                   0 to calculate the mean connectivity for all the regions:\n\
                   mat.shape = (seedRegionsNb, targetRegionsNb)\n\
                   if seed region label != 0:\n\
                   mat.shape = (seedRegionLabel vertex nb, targetRegionsNb)", true );
    app.addOption( verbose, "-verbose", "for more information.", true );
    app.addOption( normalize, "-normalize", "normalize connectivity matrix", true);
    app.addOption( connectivityTextureType, "-type",
                   "Connectivity types supported: \n\
                   (1) oneSeedRegion_to_targets (2) seedVertex_to_targets\n\
                   default = oneSeedRegion_to_targets", true);
    app.addOption( targetRegionsTexR, "-targetregionstex",
                   "target regions texture for the calculation \n\
                   of the global connectivity matrix", true );
    app.addOption( connLineTargetsDensityFileName, "-outconntargets",
                   "output connectivity targets texture file name. \n\
                   (shape (1, targetsNb))", true );
    app.addOption( connMatrixFileName, "-connmatrix",
                   "output connectivity Matrix file name",
                   true );
    app.initialize();

    //Reading inputs
    if(verbose) cout << "Reading inputs :" << endl;

    inMeshAimsR.read( inAimsMesh );
    if(verbose) cout << "cortex mesh vertex nb:"
      << int(inAimsMesh.vertex().size()) << endl;
    
    if(verbose) cout << "reading seedRegionsTex: " << flush;
    seedRegionsTexR.read( seedRegionsTex );

    if(verbose) cout << "Texture dim : " << seedRegionsTex[0].nItem() << flush;
    seedRegionsNb = textureMax(seedRegionsTex);

    if(verbose) cout << ", seedRegionsNb:" << seedRegionsNb << flush;
    if(verbose) cout << "...done." << endl;
    
    vector< size_t > * labels_ptr = labelsHistogram( seedRegionsTex, seedRegionsNb, verbose );
    vector< size_t > & labels = *labels_ptr;
    
    Matrix inConnMatrixImage;
    Connectivities * connMatrixToAllMesh_ptr = 0;
    
    vector<size_t> seedVertexIndex;
    seedVertexIndex.reserve( labels[seedRegionLabel] );
    
    if(verbose)
      cout << "labels[seedRegionLabel]:" << labels[seedRegionLabel] << endl;

    for( size_t i = 0; i < seedRegionsTex[0].nItem(); ++i )
    {
      if( seedRegionsTex[0][i] == seedRegionLabel )
        seedVertexIndex.push_back(i);  
    }

    aims::SparseOrDenseMatrix AllMeshConnMatrix;
    Reader<aims::SparseOrDenseMatrix> spmreader( inConnMatrixFileName );
    spmreader.read( AllMeshConnMatrix );

    uint cortexMeshVertexNb = uint( inAimsMesh.vertex().size() );
    connMatrixToAllMesh_ptr = 0;
    
    // sparse matrix formats
    cout << "sparse format for connectivity matrix:"
      << inConnMatrixFileName << endl;
    cout << "connectivity matrix is read as "
      << ( AllMeshConnMatrix.isDense() ? "dense" : "sparse" )
      << " matrix" << endl;
    cout << "Rows: " << size_t( AllMeshConnMatrix.getSize1() )
      <<  ", Columns: " << size_t( AllMeshConnMatrix.getSize2() ) << endl;
    connMatrixToAllMesh_ptr
      = new constel::Connectivities(labels[seedRegionLabel],
        til::SparseVector<double>(AllMeshConnMatrix.getSize2()));
    Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;

    if( AllMeshConnMatrix.getSize1() == labels[seedRegionLabel] )
    {
      if (verbose)
        cout << "input connectivity matrix correspond to seed region"
          << endl;
      for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
      {
/*        connMatrixToAllMesh[i]
          = AllMeshConnMatrix.getSparseRow
            <til::SparseVector<double> >((int32_t)i);*/
      }
    }
    else
    {
      if (verbose) cout << "input connectivity matrix correspond to all the mesh"
          << endl;
      cout << "this mode is not supported. Please use a region connectivity matrix.\n";
      return EXIT_FAILURE;
    }

    if( targetRegionsTexR.fileName() != "" )
    {
      if(verbose) cout << "reading targetRegionsTex..." << flush;
      targetRegionsTexR.read( targetRegionsTex );

      if(verbose) cout << "Texture dim : " << targetRegionsTex[0].nItem() << flush;

      targetRegionsNb = textureMax( targetRegionsTex );
      if(verbose) cout << ", targetRegionsNb: " << targetRegionsNb << flush;
      if(verbose) cout << "...done." << endl;
    }
    else
    {
      throw runtime_error( "no or wrong targetRegionsTex input...!" );
    }

    if(verbose) cout << "Label name: " << seedRegionLabel << endl;
    if(verbose) cout << "Reading inputs done." << endl << endl;
    
//     /*free the matrix*/
//     if( connMatrixToAllMesh_ptr )
//       AllMeshConnMatrix = aims::SparseOrDenseMatrix();
    
    if( connectivityTextureType == "oneSeedRegion_to_targets"
      or connectivityTextureType == "seedVertex_to_targets" )
    {
      /*
      Computing densityTexture:
      */
//       if(connMatrixToAllMesh_ptr)
//       if( AllMeshConnMatrix )
//         cout << "exist:" << endl;
//       else
//         cout << "not exist !" << endl;
//       constel::Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
//       for (uint i = 0; i < cortexMeshVertexNb; ++i)
//       {
//         for (uint j = 0; j < cortexMeshVertexNb; ++j)
//         {
//           float connMatrix_value = connMatrixToAllMesh[i][j];
//           if (connMatrix_value != 0)
//           {
//             std::cout << "(i,j): connMatrix_value:" << connMatrix_value << std::endl;
//           }
//         }
//       }
      //Computing mat [region vertex][meshVertexNb]
        // Read a regions (gyri for example) labeled texture and make a Label Histogram
        // Assumes the regions (gyri) image is labeled between 1 and regionsNb, background = 0 or lower than 1
        // Count all the labels greater than regionsNb and the labels of the background
      
      
      vector<size_t> * seedVertexIndex_ptr = &seedVertexIndex;
//       Connectivities * extractAndRegroupConnMatrix_ptr = connMatrixSeedMesh_to_targetMeshTargets_regroup( connMatrixToAllMesh_ptr, targetRegionsTex, targetRegionsNb);
      Connectivities * extractAndRegroupConnMatrix_ptr = connMatrixSeedMesh_to_targetMeshTargets_regroup( AllMeshConnMatrix, targetRegionsTex, targetRegionsNb );

      AllMeshConnMatrix = SparseOrDenseMatrix();
      SparseOrDenseMatrix *mat2 = connectivitiesToSparseOrDenseMatrix( *extractAndRegroupConnMatrix_ptr );
      AllMeshConnMatrix = *mat2;
      delete mat2;
      delete extractAndRegroupConnMatrix_ptr;
      cout << "after connMatrixSeedMesh_to_targetMeshTargets_regroup, matrix size: " << AllMeshConnMatrix.getSize1() << " x " << AllMeshConnMatrix.getSize2() << endl;
      
      if (connectivityTextureType == "oneSeedRegion_to_targets")
      {
        vector<double> * totalConnSeedRegionToTargets_ptr = connMatrixSumRows( AllMeshConnMatrix );
        double sum_i = 0;
        if (normalize)
        {
          vector<double> & line_i = *totalConnSeedRegionToTargets_ptr;
          unsigned n = line_i.size();
          for( int i = 0; i<n ; ++i )
          {
            sum_i += line_i[i];
          }
          for( size_t j = 0; j<n; ++j )
          {
            line_i[j] /= sum_i;
          }
        }
        
        vector<double>  & totalConnSeedRegionToTargets = *totalConnSeedRegionToTargets_ptr;
        TimeTexture<float> * outputTargetDensityTex_ptr = oneTargetDensityTargetsRegroupTexture( totalConnSeedRegionToTargets_ptr, targetRegionsTex, 0 );
        TimeTexture<float> & outputTargetDensityTex = *outputTargetDensityTex_ptr;
       
        if(verbose)
          cout << "Writing Density texture, shape (1, meshVertexNb):" << connTextureFileName << endl;

        Writer<TimeTexture<float> > wt( connTextureFileName );
        wt.write( outputTargetDensityTex );

        if( connLineTargetsDensityFileName != "" )
        {
          TimeTexture<float> targetsDensity_tex;
          targetsDensity_tex[0].reserve( targetRegionsNb );
          for ( size_t t = 0; t < targetRegionsNb; ++t )
          {
            targetsDensity_tex[0].push_back( totalConnSeedRegionToTargets[t] );
          }
          Writer<TimeTexture<float> > wt2( connLineTargetsDensityFileName );
          wt2.write( targetsDensity_tex );
          if(verbose)
            cout << "targets density connectivity texture written...! " << endl;
        }
        else
        {
          if(verbose)
            cout << "Don't write connectivity of targets texture, shape (1, targetsNb)." << endl;
        }
        if ( connMatrixFileName != "" )
        {
          Writer<SparseOrDenseMatrix> w( connMatrixFileName );
          w.write( AllMeshConnMatrix );
        }
        else
        {
          if(verbose)
            cout << "Don't write connectivity matrix image." << endl;
        }
        
        delete totalConnSeedRegionToTargets_ptr;
        delete outputTargetDensityTex_ptr;
        
      }
      else // (connectivityTextureType == "seedVertex_to_targets")
      {
        // TODO connMatrixNormalize( SparseOrDenseMatrix )
        if( normalize )
          connMatrixNormalize( AllMeshConnMatrix );
        cout << "seedVertex_to_targets connMatrix computed." << endl;
        if( connMatrixFileName != "" )
        {
          Writer<SparseOrDenseMatrix> w( connMatrixFileName );
          w.write( AllMeshConnMatrix );
        }
        else
          cout << "Don't write connectivity matrix image." << endl;
      }
      
      delete labels_ptr;
    }
    
    delete connMatrixToAllMesh_ptr;
    
    return EXIT_SUCCESS;
  }
  catch( carto::user_interruption & )
  {
  // Exceptions thrown by command line parser (already handled, simply exit)
  }
  catch( exception & e )
  {
    cerr << e.what() << endl;
  }
  return EXIT_FAILURE;
}
