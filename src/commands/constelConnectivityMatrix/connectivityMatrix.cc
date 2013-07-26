#include <constellation/bundleTools.h>
#include <constellation/bundleLoader.h>
#include <constellation/connMatrix.h>
#include <constellation/connMatrixTools.h>
#include <constellation/sparseMatrixSmoothing.h>
#include <constellation/connConversion.h>
#include <aims/mesh/texturetools.h>
#include <aims/getopt/getopt2.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace comist;
using namespace constel;

/*Command example:

Mode mesh_to_mesh, oneSeedRegion_to_mesh:

comistConnectivityDensityTexture \
-bundles /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/data_tests/sujet04_tracking_whiteMatter_selection_gyrus13-corps_calleux.bundles \
-mesh /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/sujet04/anatomy/meshes/sujet04_white.mesh \
-type mesh_to_mesh \
-trs /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/sujet04/referential/o200_t2_corrected_TO_nobias_anatomy.trm \
-outconntex /home/pr216633/volatile/roca/tmp/densityConnMat26_09_2008_mesh_to_mesh.tex \
-verbose \

comistConnectivityDensityTexture \
-bundles /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/data_tests/sujet04_tracking_whiteMatter_selection_gyrus13-corps_calleux.bundles \
-mesh /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/sujet04/anatomy/meshes/sujet04_white.mesh \
-type oneSeedRegion_to_mesh \
-seedlabel 13 \
-seedregionstex /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/sujet04/segment/sujet04_white_gyri2.tex \
-trs /home/pr216633/volatile/roca/thesis/diffusion/data_base_diffusion/sujet04/referential/o200_t2_corrected_TO_nobias_anatomy.trm \
-outconntex /home/pr216633/volatile/roca/tmp/densityConnMat26_09_2008_oneSeedRegion_to_mesh.tex \
-verbose \

*/


Connectivities* makeConnMatrix_intersection(
  const string & connMatrixComputingType,
  const string & bundleFilename,
  const AimsSurfaceTriangle & inAimsMesh,
  const string & motionName,
  double distthresh,
  float meshClosestPoint_maxDistance,
  AimsData<short> & roisMask,
  uint length_min,
  uint length_max,
  bool verbose )
{
  uint cortexMeshVertexNb = uint(inAimsMesh.vertex().size());
  Connectivities* connMatrixToAllMesh_ptr
    = new Connectivities(cortexMeshVertexNb,
                               til::SparseVector<double>(cortexMeshVertexNb));
  BundleInteractionReader bundleInteractionReader(bundleFilename);
  BundleProducer *finalProducer = &bundleInteractionReader;
  if ( verbose ) cout << "done" << endl;
  rc_ptr< BundleMotion > bundleMotion;
  if ( motionName.size() ) {
    if ( verbose ) cout << "Creating motion filter: " << flush;
    bundleMotion.reset( new BundleMotion( motionName ) );
    finalProducer->addBundleListener( *bundleMotion );
    finalProducer = bundleMotion.get();
    if ( verbose ) cout << "done" << endl;
  }
  MemAntBundleListener memAntBundleListener(bundleInteractionReader);
  CurvilinearAbscissaBundleListener curvilinearAbscissaBundleListener(bundleInteractionReader);

  BundleListener *meshIntersectionBundleListener = 0;
  if ( true ) // distthresh == 0.0)
  {
    if (verbose) std::cout << "nosmoothing" << std::endl;
    if ( connMatrixComputingType == "meshintersectionpointfast")
    {
      delete meshIntersectionBundleListener;
      meshIntersectionBundleListener
        = new MeshIntersectionNoSmoothingFasterBundleListener( inAimsMesh,
          roisMask, bundleInteractionReader, meshClosestPoint_maxDistance );
      static_cast<MeshIntersectionNoSmoothingFasterBundleListener *>(
        meshIntersectionBundleListener )->setMeshIdentity(0);
    }
    else
    {
      delete meshIntersectionBundleListener;
      meshIntersectionBundleListener
        = new MeshIntersectionNoSmoothingBundleListener( inAimsMesh,
          bundleInteractionReader, meshClosestPoint_maxDistance );
      static_cast<MeshIntersectionNoSmoothingBundleListener *>(
      meshIntersectionBundleListener )->setMeshIdentity(0);
    }
  }
  else // distthresh != 0
  {
    //   MeshIntersectionBundleListener meshIntersectionBundleListener(inAimsMesh, bundleInteractionReader, distthresh, meshClosestPoint_maxDistance);
    meshIntersectionBundleListener
      = new MeshIntersectionBundleListener(inAimsMesh, bundleInteractionReader,
        distthresh, meshClosestPoint_maxDistance);
    static_cast<MeshIntersectionBundleListener *>(
      meshIntersectionBundleListener )->setMeshIdentity(0);
  }
  MeshConnectionBundleListener meshConnectionBundleListener(bundleInteractionReader, 0);

  finalProducer->addBundleListener(curvilinearAbscissaBundleListener);
  finalProducer->addBundleListener(*meshIntersectionBundleListener);
  finalProducer->addBundleListener(meshConnectionBundleListener);
  finalProducer->addBundleListener(memAntBundleListener);
  if (verbose) std::cout << "Reading bundleInteractionReader..." << std::endl;
  bundleInteractionReader.read();

  //Get connections:
  ////Cortex:
  if(verbose) std::cout << "Get cortex to cortex Connections" << std::endl;
  boost::shared_ptr<BundleConnections> cortexConnections_ptr;
  boost::shared_ptr<ConnectionsLength> cortexConnectionsLength_ptr;
  cortexConnections_ptr = meshConnectionBundleListener.getBundleMeshConnections();
  cortexConnectionsLength_ptr = meshConnectionBundleListener.getBundleMeshConnectionsLength();

  //fillConnectivity Matrix:
  if(verbose) std::cout << "Fill cortex to cortex Connectivity matrix" << std::endl;
  BundleConnections & cortexConnections = *cortexConnections_ptr;
  ConnectionsLength & cortexConnectionsLength = *cortexConnectionsLength_ptr;
  if (length_min == 0 and length_max == 0)
  {
    // fillconnMatrixWithConnections(connMatrixToAllMesh_ptr, cortexConnections, wthresh, distthresh);
    fillconnMatrixWithConnections(connMatrixToAllMesh_ptr, cortexConnections, 0, 0);
  }
  else
  {
    // fillconnMatrixWithConnectionsPlusLength(connMatrixToAllMesh_ptr, cortexConnections, wthresh, distthresh, length_min, length_max, cortexConnectionsLength);
    // fillconnMatrixWithConnectionsPlusFloatLengthWeight(connMatrixToAllMesh_ptr, cortexConnections, wthresh, distthresh, length_min, length_max, cortexConnectionsLength);
    fillconnMatrixWithConnectionsPlusLengthWeight(connMatrixToAllMesh_ptr, cortexConnections, 0, 0, length_min, length_max, cortexConnectionsLength);
  }

  delete meshIntersectionBundleListener;
  return connMatrixToAllMesh_ptr;
}


Connectivities* makeConnMatrix_closestPoint(
  const string & bundleFilename,
  const AimsSurfaceTriangle & inAimsMesh, const string & motionName,
  bool verbose )
{
  BundleLoader loader;
  BundleReader bundleReader(bundleFilename);
  bundleReader.addBundleListener(loader);
  if (verbose)  std::cout << "Reading fibers..." << std::flush;
  bundleReader.read();
  rc_ptr<Fibers> pfibers;
  pfibers = loader.getFibers();
  MotionReader mreader(motionName);
  Motion motion;
  mreader.read(motion);
  // Computing complete connectivity matrix [meshVertexNb][meshVertexNb]
//       connMatrixToAllMesh_ptr = connMatrix(*pfibers, inAimsMesh, distthresh, wthresh, motion);
  Connectivities* connMatrixToAllMesh_ptr
    = connMatrix(*pfibers, inAimsMesh, 0, 0, motion);

  return connMatrixToAllMesh_ptr;
}


void makeConnectivityTexture_seedConnectionDensity(
  Connectivities* connMatrixToAllMesh_ptr,
  const string & connTextureFileName,
  TimeTexture<short> & seedRegionsTex, size_t seedRegionLabel,
  bool verbose )
{
  /*
  Computing densityTexture:
  */
  if(connMatrixToAllMesh_ptr)
  {
    if (verbose) std::cout << "exist:" << std::endl;
  }
  else
  {
    if (verbose) std::cout << "not exist !" << std::endl;
  }
  TimeTexture<float> outputTargetDensityTex = meshDensityTexture(connMatrixToAllMesh_ptr);
  if (verbose) std::cout << "Writing connectivity Density texture:" << connTextureFileName << std::endl;
  TimeTexture<float> outputTargetDensityTex_def;
  size_t cortexMeshVertexNb = outputTargetDensityTex[0].nItem();
  outputTargetDensityTex_def[0].reserve(cortexMeshVertexNb);
  for (std::size_t v = 0; v < cortexMeshVertexNb; ++v)
  {
    outputTargetDensityTex_def[0].push_back(0);
  }
  for (std::size_t v = 0; v < cortexMeshVertexNb; ++v)
  {
    if (seedRegionsTex[0][v] == seedRegionLabel)
    {
      outputTargetDensityTex_def[0][v] = outputTargetDensityTex[0][v];
    }
  }
  til::aimswrite(outputTargetDensityTex_def, connTextureFileName);
}


void makeConnectivityTexture_seedMeanConnectivityProfile(
  Connectivities* connMatrixToAllMesh_ptr,
//                                        const SparseOrDenseMatrix & connMatrixToAllMesh                  
  const string & connTextureFileName,
  const TimeTexture<short> & seedRegionsTex, size_t seedRegionLabel,
  size_t seedRegionsNb, double distthresh, double wthresh, bool logOption,
  const string & logFile, const string & connMatrixFileName,
  const AimsSurfaceTriangle & inAimsMesh, bool normalize,
  const string & connMatrixFormat, const string & seedRegionVertexIndexType,
  const string & seedRegionVertexIndexFileName,
  const vector<string> & connTextureToTargetMeshesFileNames, int meshes_nb,
  vector<Connectivities*> & connMatrixCortexToMesh_ptr_vector,
  bool verbose )
{
  if (seedRegionLabel <= 0 or seedRegionLabel > seedRegionsNb)
  {
    throw runtime_error("no or wrong seedRegionLabel...!");
  }
  //Computing mat [region vertex][meshVertexNb]
  // Read a regions (gyri for example) labeled texture and make a Label Histogram
  // Assumes the regions (gyri) image is labeled between 1 and regionsNb, background = 0 or lower than 1
  // Count all the labels greater than regionsNb and the labels of the background
  std::vector< std::size_t > * labels_ptr = labelsHistogram(seedRegionsTex, seedRegionsNb, verbose);
  std::vector< std::size_t > & labels = *labels_ptr;

  vector<size_t> * seedVertexIndex;
  std::size_t seedRegionLabelVertexNb = labels[seedRegionLabel];
  Connectivities * extractConnMatrix_ptr = connMatrixReducedFromRegion( connMatrixToAllMesh_ptr, seedRegionsTex, seedRegionLabel, seedRegionLabelVertexNb, & seedVertexIndex );
  if( distthresh != 0 )
    sparseMatrixDiffusionSmoothing( extractConnMatrix_ptr, inAimsMesh,
      wthresh, distthresh, seedRegionsTex, seedRegionLabel );

  Connectivities & extractConnMatrix = *extractConnMatrix_ptr;
//   AllMeshConnMatrix = SparseOrDenseMatrix();
  SparseOrDenseMatrix *mat = connectivitiesToSparseOrDenseMatrix( *extractConnMatrix_ptr );
//   AllMeshConnMatrix = *mat;
  if(normalize)
    connMatrixNormalize( *mat );
  if(connMatrixFileName != "" )
  {
    Writer<SparseOrDenseMatrix> w( connMatrixFileName );
    w.write( *mat );
    if (logOption and logFile != "")
    {
      if (verbose) std::cout << "Computing ln(1+matrix) and storing resulting file in " << connMatrixFileName << "..." << std::endl;
      std::size_t colNb = extractConnMatrix[0].size();
      std::size_t rowsNb = extractConnMatrix.size();
      Connectivities::iterator il, el = extractConnMatrix.end();
      Connectivity::sparse_iterator ic, ec;
      for( il=extractConnMatrix.begin(); il!=el; ++il )
      {
        for( ic=il->sparse_begin(), ec=il->sparse_end(); ic!=ec; ++ic )
        {
          float connval = ic->second;
          if (connval > 0)
          {
            float new_connval = log(1 + connval);
            ic->second = new_connval;
          }
        }
      }

      if (connMatrixFormat == "binar_sparse" or connMatrixFormat == "ascii")
      {
        if(verbose) std::cout << "Writing log seed region connMatrix ima:" << connMatrixFileName << std::endl;
        writeConnectivities( *extractConnMatrix_ptr, logFile,
                             connMatrixFormat=="ascii" );
      }
      else//saving connectivity matrix in .ima (aims fmt)
      {
        writeAimsFmtConnMatrix(extractConnMatrix_ptr,logFile);
      }
    }
    //       Write region vertex indexes (ascii format only)
    std::size_t seedVertexIndex_size = (*seedVertexIndex).size();
    if ((*seedVertexIndex).size()!=0)
    {
      if (seedRegionVertexIndexType == "text" or seedRegionVertexIndexType == "both")
      {
        std::fstream findex;
        ostringstream s;
        s << seedRegionVertexIndexFileName << ".txt";
        findex.open(s.str().c_str(), std::fstream::out);
        for (std::size_t i = 0; i < labels[seedRegionLabel]; ++i)
        {
          findex << (*seedVertexIndex)[i] << std::endl;
        }
      }
      else if (seedRegionVertexIndexType == "texture" or seedRegionVertexIndexType == "both")
      {
        TimeTexture< unsigned int > seedRegionVertexIndexTex;
        if (verbose) std::cout << "seedVertexIndex_size:" << seedVertexIndex_size << std::endl;
        seedRegionVertexIndexTex[0].reserve(seedVertexIndex_size);
        for (std::size_t vertex = 0; vertex < seedVertexIndex_size; vertex++)
        {
          unsigned int val = (*seedVertexIndex)[vertex];
          seedRegionVertexIndexTex[0].push_back((*seedVertexIndex)[vertex]);
        }
        if (verbose) std::cout << "seedRegionVertexIndexTex0.nItem():" << seedRegionVertexIndexTex.nItem() << ", " << seedRegionVertexIndexTex[0].nItem() << std::endl;
        ostringstream s;
        s << seedRegionVertexIndexFileName << ".tex";
        Writer<TimeTexture<unsigned int > > w(s.str());
        w.write(seedRegionVertexIndexTex);
      }
      if (verbose) std::cout << "End Writing seedVertexIndex in " << seedRegionVertexIndexType << " format." << std::endl;
    }
    else
    {
      if (verbose) std::cout << "seedVertexIndex size = 0" << std::endl;
    }
  }
  /*
  Computing densityTexture:
    outputDensityTex: texture of the connection density of the seed region towards all the other vertex of the mesh
    outputTargetDensityTex: texture of the connection density of the entire mesh towards the seed region
  */
//         TimeTexture<float> outputDensityTex = densityTexture(extractConnMatrix_ptr, *seedVertexIndex);
//         ostringstream s;
//         s << outputDirname << "DensityConnection_region" << seedRegionLabel << "_" << outputString;
//         til::aimswrite(outputDensityTex, s.str());
  TimeTexture<float> outputTargetDensityTex = meshDensityTexture(extractConnMatrix_ptr);
  if (verbose) std::cout << "Writing mean connectivity profile texture:" << connTextureFileName << std::endl;
  til::aimswrite(outputTargetDensityTex, connTextureFileName);
  delete extractConnMatrix_ptr;
  delete seedVertexIndex;
  int connTextureToTargetMeshesFileNames_size = connTextureToTargetMeshesFileNames.size();
  if( connTextureToTargetMeshesFileNames_size == meshes_nb )
  {
    for( int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel )
    {
      if( connTextureToTargetMeshesFileNames[meshLabel] != "" )
      {
        vector<size_t> * seedVertexIndex;
        SparseOrDenseMatrix *mat = connectivitiesToSparseOrDenseMatrix( *connMatrixCortexToMesh_ptr_vector[meshLabel] );
        Connectivities * extractConnMatrix_ptr = connMatrixRegionExtract( *mat, seedRegionsTex, seedRegionLabel, seedRegionLabelVertexNb, & seedVertexIndex );
        SparseOrDenseMatrix *mat2 = connectivitiesToSparseOrDenseMatrix( *extractConnMatrix_ptr );
        if(normalize)
          connMatrixNormalize( *mat2 );
        TimeTexture<float> outputTargetDensityTex = meshDensityTexture( *mat2 );
        if(verbose) std::cout << "Writing mean connectivity profile texture:" << connTextureToTargetMeshesFileNames[meshLabel] << std::endl;
        til::aimswrite(outputTargetDensityTex, connTextureToTargetMeshesFileNames[meshLabel]);
        delete extractConnMatrix_ptr;
        delete seedVertexIndex;
      }
    }
  }

  delete labels_ptr;
}



int main( int argc, char* argv[] )
{
  try
  {
    std::string bundleFilename;
    Reader<AimsSurfaceTriangle> inMeshAimsR;
    std::string connMatrixComputingType = "meshclosestpoint";
    Reader< TimeTexture<short> > seedRegionsTexR;
    AimsSurfaceTriangle inAimsMesh;
    std::vector< std::string > inTargetMeshesAimsR_files;
    TimeTexture<short> seedRegionsTex;
    std::size_t seedRegionLabel = 0;
    std::size_t seedRegionsNb = 0;
    double distthresh = 5.0;
    double wthresh = 1.0;
    std::string motionName = "";
    std::string connTextureFileName;
    std::vector<std::string> connTextureToTargetMeshesFileNames;
    connTextureToTargetMeshesFileNames.reserve(1);
    connTextureToTargetMeshesFileNames.push_back("");
    std::string connMatrixFileName = "";
    std::string connMatrixFormat = "binar_sparse";
    bool verbose = false;
    bool normalize = false;
    double meshesDistanceThreshold = 1.0;
    float meshClosestPoint_maxDistance = 5.0;//in mm
    std::string connectivityTextureType = "seed_connection_density";
    std::string seedRegionVertexIndexFileName;
    std::string  seedRegionVertexIndexType = "";
    uint length_min = 0;
    uint length_max = 0;
    bool logOption = false;
    std::string logFile = "";
    Reader<AimsData<short> > roisMaskR;
    AimsData<short> roisMask;
    
    AimsApplication app( argc, aims_const_hack(argv),
                         "computes and writes connection density texture(s).\n\
                         Modes supported: (1) mesh_to_mesh (2) oneSeedRegion_to_mesh" );
    app.addOption( bundleFilename, "-bundles", "input bundles" );
    app.addOption( inMeshAimsR, "-mesh", "input mesh" );
    app.addOptionSeries( inTargetMeshesAimsR_files, "-targets",
                         "input target meshes" );
    app.addOption( connTextureFileName, "-outconntex",
                   "output mean connectivity texture file name" );
    app.addOptionSeries( connTextureToTargetMeshesFileNames, "-outtargetstex",
                         "output mean connectivity texture file name\n\
                         (to targets meshes)");
    app.addOption( motionName, "-trs", "transform from t2 to anat" );
    app.addOption( connMatrixComputingType, "-matrixcompute",
                   "type of computing the connectivity matrix supported:\n\
                   (1) meshclosestpoint (2) meshintersectionpoint (3) meshintersectionpointfast\n\
                   default = meshclosestpoint", true );
    app.addOption( distthresh, "-dist",
                   "dist for the neighborhood around each vertex for smoothing the connectivity matrix\n\
                   default = 5.0", true );
    app.addOption( wthresh, "-wthresh",
                   "weight threshold for thresholding the connectivity matrix, default = 1.0", true );
    app.addOption( meshClosestPoint_maxDistance, "-distmax",
                   "mesh closest point minimum distance, default = 25.", true );
    app.addOption( seedRegionsTexR, "-seedregionstex",
                   "input region texture. default = all the mesh", true );
    app.addOption( seedRegionLabel, "-seedlabel",
                   "input seed region label. 0 to calculate the mean connectivity for all the regions : connMatrix.shape = (seedRegionsNb, targetRegionsNb); = (seedRegionLabel vertex nb, targetRegionsNb) if seed region label != 0", true );
    app.addOption( verbose, "-verbose", "show as much information as possible", true );
    app.addOption( normalize, "-normalize", "normalize connectivity matrix", true);
    app.addOption( connectivityTextureType, "-type",
                   "connectivity type: seed_mean_connectivity_profile, default = seed_connection_density");
    app.addOption( connMatrixFormat, "-connfmt",
                   "input conn matrix format, .ima or binar_sparse; default = binar_sparse",true );
    app.addOption( connMatrixFileName, "-connmatrix",
                   "connectivity matrix filename of the seed region, \n\
                   size = (seedRegionVertexNb,meshVertexNb)",true);
    app.addOption( seedRegionVertexIndexType, "-vertexindextype",
                   "format for seedRegionVertexIndexFile:\n\
                   =texture(.tex) or =text(.txt) or =both(.tex and .txt), default = "" ", true );
    app.addOption( seedRegionVertexIndexFileName, "-seedvertexindex",
                   "seedlabel region vertex indexes file name", true );
    app.addOption( length_min, "-lmin", "length min for fibers", true);
    app.addOption( length_max, "-lmax", "length max for fibers", true);
    app.addOption( logOption, "-log", "log option, if = True,\n\
                   applies the ln(1 + x) to the connectivity matrix, default = False", true);
    app.addOption( logFile, "-logfile",
                   "log filename, in case of logOption = 1, store log matrix in another file.\n\
                   by default, store log matrix in connMatrixFileName", true);
    app.addOption( roisMaskR, "-roimask",
                   "input region roi mask: ribbon around the mesh, in the case of  connMatrixComputingType = meshintersectionpointfast", true );
    app.initialize();

    //Reading inputs
    if(verbose) std::cout << "Reading inputs :" << std::endl;
    if(verbose) std::cout << "input bundles : " << bundleFilename << std::endl;
    if(verbose) std::cout << "dist : " << distthresh << std::endl;
    if(verbose) std::cout << "wthresh : " << wthresh << std::endl;

    inMeshAimsR.read( inAimsMesh );
    if(verbose) std::cout << "cortex mesh vertex nb:" << int(inAimsMesh.vertex().size()) << std::endl;
    
    uint cortexMeshVertexNb = uint(inAimsMesh.vertex().size());
    int meshes_nb = inTargetMeshesAimsR_files.size();
    std::vector<AimsSurfaceTriangle> inAimsTargetMeshes(meshes_nb);
    
    //Reading target meshes:
    if(inTargetMeshesAimsR_files.empty() )
    {
      if(verbose) std::cout << "No target meshes" << std::endl;
    }
    else
    {
      if (verbose) std::cout << "Reading target meshes:" << std::endl;
//       meshes_nb = inTargetMeshesAimsR_files.size();
      if(verbose) std::cout << "meshes_nb:" << meshes_nb << std::endl;
      if (verbose) std::cout << "inAimsTargetMeshes.size(): " << inAimsTargetMeshes.size() << std::endl;
      for (int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel)
      {
        if(verbose) std::cout << meshLabel << ":" << inTargetMeshesAimsR_files[meshLabel] << std::endl;
        Reader<AimsSurfaceTriangle> inAimsTargetMeshR(inTargetMeshesAimsR_files[meshLabel]);
        inAimsTargetMeshR.read( inAimsTargetMeshes[meshLabel]);
        if(verbose) std::cout << "ok..." << std::flush;
      }
      if (verbose) std::cout << "inAimsTargetMeshes.size(): " << inAimsTargetMeshes.size() << std::endl;
      if (verbose) std::cout << "inAimsTargetMesh0.vertex().size():" << inAimsTargetMeshes[0].vertex().size() << std::endl;
      if(verbose) std::cout << "Done." << std::endl;
    }
    
    if(seedRegionsTexR.fileName()!= "")
    {
      if(verbose) std::cout << "reading seedRegionsTex: " << flush;
      seedRegionsTexR.read( seedRegionsTex );
      if(verbose) std::cout << "Texture dim : " << seedRegionsTex[0].nItem() << std::flush;
      seedRegionsNb = textureMax(seedRegionsTex);
      if(verbose) std::cout << ", seedRegionsNb:" << seedRegionsNb << std::flush;
      if(verbose) std::cout << "...done." << std::endl;
    }
    else
    {
      throw runtime_error("no or wrong seedRegionsTex input...!");
    }
    if(roisMaskR.fileName()!= "")
    {
      if(verbose) std::cout << "reading roi mask Volume: " << flush;
      roisMaskR.read( roisMask );
      if(verbose) std::cout << "...done." << std::endl;
    }
    
    if (verbose)
    {
      std::cout << "seedRegionLabel : " << seedRegionLabel << std::endl;
    }
    if(verbose) std::cout << "Reading inputs done." << std::endl << std::endl;
    Connectivities * connMatrixToAllMesh_ptr = 0;
    // loading fibers
    if(verbose) std::cout << "Loading fibers..." << std::flush;

    std::vector<Connectivities* > connMatrixCortexToMesh_ptr_vector;
    connMatrixCortexToMesh_ptr_vector.reserve(meshes_nb);
    if (connMatrixComputingType == "meshintersectionpoint" ||
        connMatrixComputingType == "meshintersectionpointfast" )
      connMatrixToAllMesh_ptr = makeConnMatrix_intersection(
        connMatrixComputingType, bundleFilename, inAimsMesh, motionName,
        distthresh, meshClosestPoint_maxDistance, roisMask, length_min,
        length_max, verbose );
    else //connMatrixComputingType == "meshclosestpoint"
      connMatrixToAllMesh_ptr = makeConnMatrix_closestPoint(
        bundleFilename, inAimsMesh, motionName, verbose );

    if(verbose) std::cout << "OK" << std::endl;

    
    
    if (connectivityTextureType == "seed_connection_density")
      makeConnectivityTexture_seedConnectionDensity( connMatrixToAllMesh_ptr,
        connTextureFileName, seedRegionsTex, seedRegionLabel, verbose );
    else if (connectivityTextureType == "seed_mean_connectivity_profile")
      makeConnectivityTexture_seedMeanConnectivityProfile(
        connMatrixToAllMesh_ptr, connTextureFileName, seedRegionsTex,
        seedRegionLabel, seedRegionsNb, distthresh, wthresh, logOption, logFile, connMatrixFileName, inAimsMesh, normalize, connMatrixFormat, seedRegionVertexIndexType, seedRegionVertexIndexFileName, connTextureToTargetMeshesFileNames, meshes_nb, connMatrixCortexToMesh_ptr_vector, verbose );
    delete connMatrixToAllMesh_ptr;
    for (int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel)
    {
//       delete connMatrixCortexToMesh_ptr_vector[meshLabel];
    }
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
}

