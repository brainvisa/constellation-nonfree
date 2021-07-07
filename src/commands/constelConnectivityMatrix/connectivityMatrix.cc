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
using namespace constel;

//-----------------Connectivity Matrix: intersection with mesh-----------------

Connectivities* makeConnMatrix_intersection(
    const string & connMatrixComputingType, const string & bundleFilename,
    const AimsSurfaceTriangle & inAimsMesh, const string & motionName,
    double distthresh, float meshClosestPoint_maxDistance,
    AimsData<short> & roisMask, uint length_min, uint length_max,
    bool verbose) {
  uint cortexMeshVertexNb = uint(inAimsMesh.vertex().size());
  if (verbose) cout << "Mesh vertices number: " << cortexMeshVertexNb << endl;

  Connectivities* connMatrixToAllMesh_ptr = new Connectivities(
      cortexMeshVertexNb, Connectivity(cortexMeshVertexNb));
  
  if (verbose) cout << "Loading fibers." << endl;
  BundleInteractionReader bundleInteractionReader(bundleFilename);
  BundleProducer *finalProducer = &bundleInteractionReader;
  
  rc_ptr<BundleMotion> bundleMotion;
  if (motionName.size()) {
    if (verbose) cout << "Creating motion filter." << endl;
    bundleMotion.reset(new BundleMotion(motionName));
    finalProducer->addBundleListener(*bundleMotion);
    finalProducer = bundleMotion.get();
  }

  MemAntBundleListener memAntBundleListener(bundleInteractionReader);
  CurvilinearAbscissaBundleListener curvilinearAbscissaBundleListener(
      bundleInteractionReader);

  BundleListener *meshIntersectionBundleListener = 0;
  if (true) { // no smoothing (distthresh == 0.0)
    if (verbose) cout << "No smoothing." << endl;
    if (connMatrixComputingType == "meshintersectionpointfast") {
      delete meshIntersectionBundleListener;
      meshIntersectionBundleListener
        = new MeshIntersectionNoSmoothingFasterBundleListener(
            inAimsMesh, roisMask, bundleInteractionReader,
            meshClosestPoint_maxDistance);
      static_cast<MeshIntersectionNoSmoothingFasterBundleListener *>(
          meshIntersectionBundleListener)->setMeshIdentity(0);
    } else {
      delete meshIntersectionBundleListener;
      meshIntersectionBundleListener
        = new MeshIntersectionNoSmoothingBundleListener(
            inAimsMesh, bundleInteractionReader, meshClosestPoint_maxDistance);
      static_cast<MeshIntersectionNoSmoothingBundleListener *>(
          meshIntersectionBundleListener)->setMeshIdentity(0);
    }
  } else { // smoothing (distthresh != 0.0)
    meshIntersectionBundleListener = new MeshIntersectionBundleListener(
        inAimsMesh, bundleInteractionReader, distthresh,
        meshClosestPoint_maxDistance);
    static_cast<MeshIntersectionBundleListener *>(
        meshIntersectionBundleListener)->setMeshIdentity(0);
  }

  MeshConnectionBundleListener meshConnectionBundleListener(
      bundleInteractionReader, 0);

  finalProducer->addBundleListener(curvilinearAbscissaBundleListener);
  finalProducer->addBundleListener(*meshIntersectionBundleListener);
  finalProducer->addBundleListener(meshConnectionBundleListener);
  finalProducer->addBundleListener(memAntBundleListener);

  if (verbose) cout << "Reading bundleInteractionReader." << endl;
  bundleInteractionReader.read();

  //Get connections:
  ////Cortex:
  if (verbose) cout << "- Get cortex to cortex Connections." << endl;

  boost::shared_ptr<BundleConnections> cortexConnections_ptr;
  boost::shared_ptr<ConnectionsLength> cortexConnectionsLength_ptr;

  cortexConnections_ptr
    = meshConnectionBundleListener.getBundleMeshConnections();

  cortexConnectionsLength_ptr
    = meshConnectionBundleListener.getBundleMeshConnectionsLength();

  //fillConnectivity Matrix:
  if (verbose) cout << "- Fill cortex to cortex Connectivity matrix" << endl;

  BundleConnections & cortexConnections = *cortexConnections_ptr;
  ConnectionsLength & cortexConnectionsLength = *cortexConnectionsLength_ptr;

  if (length_min == 0 and length_max == 0)
    fillconnMatrixWithConnections(connMatrixToAllMesh_ptr,
                                  cortexConnections, 0, 0);
  else
    fillconnMatrixWithConnectionsPlusLengthWeight(
        connMatrixToAllMesh_ptr, cortexConnections, 0, 0, length_min,
        length_max, cortexConnectionsLength);

  delete meshIntersectionBundleListener;
  return connMatrixToAllMesh_ptr;
}

//-----------------Connectivity Matrix: closest point with mesh----------------

Connectivities* makeConnMatrix_closestPoint(
    const string & bundleFilename,
    const AimsSurfaceTriangle & inAimsMesh, const string & motionName,
    bool verbose) {
  BundleLoader loader;
  BundleReader bundleReader(bundleFilename);
  bundleReader.addBundleListener(loader);

  if (verbose) cout << "Reading fibers..." << flush;

  bundleReader.read();
  rc_ptr<Fibers> pfibers;
  pfibers = loader.getFibers();
  MotionReader mreader(motionName);
  Motion motion;
  mreader.read(motion);

  // Computing complete connectivity matrix [meshVertexNb][meshVertexNb]
  Connectivities* connMatrixToAllMesh_ptr = connMatrix(
      *pfibers, inAimsMesh, 0, 0, motion);

  return connMatrixToAllMesh_ptr;
}

//-----------------Texture: Seed Connection Density----------------------------

void makeConnectivityTexture_seedConnectionDensity(
    Connectivities* connMatrixToAllMesh_ptr,
    const string & connTextureFileName,
    TimeTexture<short> & seedRegionsTex, int seedRegionLabel,
    bool verbose) {
  // Computing densityTexture:
  if (connMatrixToAllMesh_ptr)
  {
    if (verbose) cout << "exist:" << endl;
  }
  else
  {
    if (verbose) cout << "not exist !" << endl;
  }

  TimeTexture<float> outputTargetDensityTex = meshDensityTexture(
      connMatrixToAllMesh_ptr);

  if (verbose)
    cout << "Writing connectivity Density texture:" << connTextureFileName <<
    endl;

  TimeTexture<float> outputTargetDensityTex_def;
  size_t cortexMeshVertexNb = outputTargetDensityTex[0].nItem();
  outputTargetDensityTex_def[0].reserve(cortexMeshVertexNb);

  for (size_t v = 0; v < cortexMeshVertexNb; ++v) {
    outputTargetDensityTex_def[0].push_back(0);
  }

  for (size_t v = 0; v < cortexMeshVertexNb; ++v) {
    if (seedRegionsTex[0][v] == seedRegionLabel) {
      outputTargetDensityTex_def[0][v] = outputTargetDensityTex[0][v];
    }
  }

  Writer<TimeTexture<float> > w( connTextureFileName );
  w.write( outputTargetDensityTex_def );
}

//-----------------Texture: Seed Mean Connectivity Profile---------------------

void makeConnectivityTexture_seedMeanConnectivityProfile(
    Connectivities* connMatrixToAllMesh_ptr,
    // const SparseOrDenseMatrix & connMatrixToAllMesh                  
    const string & connTextureFileName,
    const TimeTexture<short> & seedRegionsTex,
    int seedRegionLabel,
    int maxLabel,
    double distthresh, double wthresh,
    bool logOption, const string & logFile,
    const string & connMatrixFileName,
    const AimsSurfaceTriangle & inAimsMesh,
    bool normalize, const string & connMatrixFormat,
    const string & seedRegionVertexIndexType,
    const string & seedRegionVertexIndexFileName,
    const vector<string> & connTextureToTargetMeshesFileNames,
    int meshes_nb,
    vector<Connectivities*> & connMatrixCortexToMesh_ptr_vector,
    bool verbose)
{

  /*
  Computing mat [region vertex][meshVertexNb]
  Read a regions (gyri for example) labeled texture and make a Label Histogram.
  */
  int minLabel = textureMin(seedRegionsTex);
  if (seedRegionLabel < minLabel or seedRegionLabel > maxLabel)
    throw runtime_error("No or wrong seedRegionLabel.");
  map<short, size_t> *labels_ptr = labelsHistogram(seedRegionsTex,
                                                    maxLabel,
                                                    minLabel,
                                                    verbose);
  map<short, size_t> &labels = *labels_ptr;

  vector<size_t> *seedVertexIndex;
  int seedRegionLabelVertexNb = labels[seedRegionLabel];

  rc_ptr<Connectivities> extractConnMatrix_ptr(
    connMatrixReducedFromRegion(connMatrixToAllMesh_ptr, seedRegionsTex, 
                                seedRegionLabel, seedRegionLabelVertexNb, 
                                & seedVertexIndex));
  
  if (distthresh != 0)
    sparseMatrixDiffusionSmoothing(extractConnMatrix_ptr, inAimsMesh,
                                   wthresh, distthresh,
                                   seedRegionsTex, seedRegionLabel);

  Connectivities &extractConnMatrix = *extractConnMatrix_ptr;
  //  AllMeshConnMatrix = SparseOrDenseMatrix();
  SparseOrDenseMatrix *mat = connectivitiesToSparseOrDenseMatrix(
      *extractConnMatrix_ptr);
  //  AllMeshConnMatrix = *mat;
  
  if (normalize)
    connMatrixNormalize(*mat);
  
  if (connMatrixFileName != "")
  {
    Writer<SparseOrDenseMatrix> w(connMatrixFileName);
    w.write(*mat);
    if (logOption and logFile != "")
    {
      if (verbose)
        cout << "Computing ln(1+matrix) and storing resulting file in " <<
        connMatrixFileName << "..." << endl;
      
      //size_t colNb = extractConnMatrix[0].size();
      //size_t rowsNb = extractConnMatrix.size();
      
      Connectivities::iterator il, el = extractConnMatrix.end();
      Connectivity::iterator ic, ec;
      
      for (il=extractConnMatrix.begin(); il!=el; ++il)
      {
        for (ic=il->begin(), ec=il->end(); ic!=ec; ++ic)
        {
          float connval = ic->second;
          
          if (connval > 0)
          {
            float new_connval = log(1 + connval);
            ic->second = new_connval;
          }
        }
      }

      if (connMatrixFormat == "binar_sparse" || connMatrixFormat == "ascii")
      {
        if (verbose)
          cout << "Writing log seed region connMatrix ima:" <<
          connMatrixFileName << endl;
        
        writeConnectivities(*extractConnMatrix_ptr, logFile,
                            connMatrixFormat=="ascii");
      }
      else
      {  // saving connectivity matrix in .ima (aims fmt)
        writeAimsFmtConnMatrix(extractConnMatrix_ptr.get(),logFile);
      }
    }
    
    //  Write region vertex indexes (ascii format only)
    size_t seedVertexIndex_size = (*seedVertexIndex).size();
    if ((*seedVertexIndex).size() != 0)
    {
      if (seedRegionVertexIndexType == "text"
          || seedRegionVertexIndexType == "both")
      {
        fstream findex;
        ostringstream s;
        s << seedRegionVertexIndexFileName << ".txt";
        findex.open(s.str().c_str(), fstream::out);
        
        for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
        {
          findex << (*seedVertexIndex)[i] << endl;
        }
      }
      else if( seedRegionVertexIndexType == "texture"
               || seedRegionVertexIndexType == "both" )
      {
        TimeTexture< unsigned int > seedRegionVertexIndexTex;
        
        if (verbose)
          cout << "seedVertexIndex_size:" << seedVertexIndex_size << endl;
        
        seedRegionVertexIndexTex[0].reserve(seedVertexIndex_size);
        
        for (size_t vertex = 0; vertex < seedVertexIndex_size; vertex++)
        {
          seedRegionVertexIndexTex[0].push_back((*seedVertexIndex)[vertex]);
        }

        if (verbose)
          cout << "seedRegionVertexIndexTex0.nItem():" <<
          seedRegionVertexIndexTex.nItem() << ", " <<
          seedRegionVertexIndexTex[0].nItem() << endl;
        
        ostringstream s;
        s << seedRegionVertexIndexFileName << ".tex";
        Writer<TimeTexture<unsigned int > > w(s.str());
        w.write(seedRegionVertexIndexTex);
      }
      if (verbose)
        cout << "End Writing seedVertexIndex in " <<
        seedRegionVertexIndexType << " format." << endl;
    }
    else
    {
      if (verbose) cout << "seedVertexIndex size = 0" << endl;
    }
  }

  /*
  Computing densityTexture:
  - outputDensityTex: texture of the connection density of the seed region
    towards all the other vertex of the mesh
  - outputTargetDensityTex: texture of the connection density of the entire
    mesh towards the seed region
  */
  TimeTexture<float> outputTargetDensityTex = meshDensityTexture(
    extractConnMatrix_ptr.get());

  if (verbose)
    cout << "Writing mean connectivity profile texture:" <<
    connTextureFileName << endl;

  Writer<TimeTexture<float> > w( connTextureFileName );
  w.write( outputTargetDensityTex );
  extractConnMatrix_ptr.reset(); // delete

  delete seedVertexIndex;

  int connTextureToTargetMeshesFileNames_size
    = connTextureToTargetMeshesFileNames.size();

  if (connTextureToTargetMeshesFileNames_size == meshes_nb)
  {
    for (int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel)
    {
      if (connTextureToTargetMeshesFileNames[meshLabel] != "")
      {
        vector<size_t> * seedVertexIndex;

        SparseOrDenseMatrix *mat = connectivitiesToSparseOrDenseMatrix(
            *connMatrixCortexToMesh_ptr_vector[meshLabel]);

        Connectivities * extractConnMatrix_ptr = connMatrixRegionExtract(
            *mat, seedRegionsTex, seedRegionLabel, seedRegionLabelVertexNb,
            & seedVertexIndex);

        SparseOrDenseMatrix *mat2 = connectivitiesToSparseOrDenseMatrix(
            *extractConnMatrix_ptr);

        if (normalize) connMatrixNormalize(*mat2);

        TimeTexture<float> outputTargetDensityTex = meshDensityTexture(*mat2);

        if (verbose)
          cout << "Writing mean connectivity profile texture:" <<
          connTextureToTargetMeshesFileNames[meshLabel] << endl;

        write( outputTargetDensityTex,
               connTextureToTargetMeshesFileNames[meshLabel] );

        delete extractConnMatrix_ptr;
        delete seedVertexIndex;
      }
    }
  }
  delete labels_ptr;
}


int main(int argc, const char* argv[])
{
  try
  {
    string bundleFilename;
    Reader<AimsSurfaceTriangle> inMeshAimsR;
    string connMatrixComputingType = "meshclosestpoint";
    Reader<TimeTexture<short> > seedRegionsTexR;
    AimsSurfaceTriangle inAimsMesh;
    vector<string> inTargetMeshesAimsR_files;
    TimeTexture<short> seedRegionsTex;
    int seedRegionLabel = 0;
    size_t maxLabel = 0;
    double distthresh = 5.0;
    double wthresh = 1.0;
    string motionName = "";
    string connTextureFileName;
    vector<string> connTextureToTargetMeshesFileNames;
    connTextureToTargetMeshesFileNames.reserve(1);
    connTextureToTargetMeshesFileNames.push_back("");
    string connMatrixFileName = "";
    string connMatrixFormat = "binar_sparse";
    bool verbose = false;
    bool normalize = false;
    //double meshesDistanceThreshold = 1.0;
    float meshClosestPoint_maxDistance = 5.0;  //in mm
    string connectivityTextureType = "seed_connection_density";
    string seedRegionVertexIndexFileName;
    string  seedRegionVertexIndexType = "";
    uint length_min = 0;
    uint length_max = 0;
    bool logOption = false;
    string logFile = "";
    Reader<AimsData<short> > roisMaskR;
    AimsData<short> roisMask;

    AimsApplication app(
        argc, argv,
        "computes and writes connection density texture(s).\n\
        Modes supported: (1) mesh_to_mesh (2) oneSeedRegion_to_mesh");

    app.addOption(
        bundleFilename, "-bundles",
        "input bundles");
    app.addOption(
        inMeshAimsR, "-mesh",
        "input mesh");
    app.addOptionSeries(
        inTargetMeshesAimsR_files, "-targets",
        "input target meshes");
    app.addOption(
        connTextureFileName, "-outconntex",
        "output mean connectivity texture file name");
    app.addOptionSeries(
        connTextureToTargetMeshesFileNames, "-outtargetstex",
        "output mean connectivity texture file name (to targets meshes)");
    app.addOption(
        motionName, "-trs",
        "transform from t2 to anat");
    app.addOption(
        connMatrixComputingType, "-matrixcompute",
        "type of computing the connectivity matrix supported: \n\
        (1) meshclosestpoint (2) meshintersectionpoint \n\
        (3) meshintersectionpointfast, default = meshclosestpoint", true);
    app.addOption(
        distthresh, "-dist",
        "dist for the neighborhood around each vertex for smoothing the \n\
        connectivity matrix, default = 5.0", true);
    app.addOption(
        wthresh, "-wthresh",
        "weight threshold for thresholding the connectivity matrix, \n\
        default = 1.0", true);
    app.addOption(
        meshClosestPoint_maxDistance, "-distmax",
        "mesh closest point minimum distance, default = 5.", true);
    app.addOption(
        seedRegionsTexR, "-seedregionstex",
        "input region texture. default = all the mesh", true);
    app.addOption(
        seedRegionLabel, "-seedlabel",
        "input seed region label. 0 to calculate the mean connectivity for \n\
        all the regions : \n\
        connMatrix.shape=(seedRegionsNb, targetRegionsNb); \n\
        = (seedRegionLabel vertex nb, targetRegionsNb) \n\
        if seed region label != 0", true);
    app.addOption(
        verbose, "-verbose",
        "show as much information as possible", true);
    app.addOption(
        normalize, "-normalize",
        "normalize connectivity matrix", true);
    app.addOption(
        connectivityTextureType, "-type",
        "connectivity type: seed_mean_connectivity_profile, \n\
        default = seed_connection_density");
    app.addOption(
        connMatrixFormat, "-connfmt",
        "input conn matrix format, .ima or binar_sparse; \n\
        default = binar_sparse",true);
    app.addOption(
        connMatrixFileName, "-connmatrix",
        "connectivity matrix filename of the seed region, \n\
        size = (seedRegionVertexNb,meshVertexNb)",true);
    app.addOption(
        seedRegionVertexIndexType, "-vertexindextype",
        "format for seedRegionVertexIndexFile: \n\
        =texture(.tex) or =text(.txt) or =both(.tex and .txt), \n\
        default = "" ", true);
    app.addOption(
        seedRegionVertexIndexFileName, "-seedvertexindex",
        "seedlabel region vertex indexes file name", true);
    app.addOption(
        length_min, "-lmin",
        "length min for fibers", true);
    app.addOption(
        length_max, "-lmax",
        "length max for fibers", true);
    app.addOption(
        logOption, "-log",
        "log option, if = True, applies the ln(1 + x) to the connectivity \n\
        matrix, default = False", true);
    app.addOption(
        logFile, "-logfile",
        "log filename, in case of logOption = 1, store log matrix in \n\
        another file. \n\
        By default, store log matrix in connMatrixFileName", true);
    app.addOption(
        roisMaskR, "-roimask",
        "input region roi mask: ribbon around the mesh, in the case of \n\
        connMatrixComputingType = meshintersectionpointfast", true);
    app.initialize();

    if (verbose)
    {
      cout << "Cortical region: " << seedRegionLabel << endl;
    }

    // Read the mesh file
    inMeshAimsR.read(inAimsMesh);

    // Define the number of target mesh files
    int meshes_nb = inTargetMeshesAimsR_files.size();

    // Target meshes (useful?)
    if (inTargetMeshesAimsR_files.empty())
    {
      if (verbose) cout << "No target meshes." << endl;
    }
    else
    {
      vector<AimsSurfaceTriangle> inAimsTargetMeshes(meshes_nb);
      if (verbose)
      {
        cout << "Meshes number: " << meshes_nb << endl;
      }
      for (int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel)
      {
        // Read the target mesh
        Reader<AimsSurfaceTriangle> inAimsTargetMeshR(
            inTargetMeshesAimsR_files[meshLabel]);
        inAimsTargetMeshR.read(inAimsTargetMeshes[meshLabel]);

        if (verbose)
        {
          cout << "Mesh: " << inTargetMeshesAimsR_files[meshLabel] << flush;
          cout << ", size of mesh: " <<
          inAimsTargetMeshes[meshLabel].vertex().size() << flush;
          cout << " ...done." << endl;
        }
      }
    }

    // Read of the starting cortical parcellation (if it exists) and give the
    // number of different regions in this parcellation
    if (seedRegionsTexR.fileName() != "")
    {
      seedRegionsTexR.read(seedRegionsTex);
      maxLabel = textureMax(seedRegionsTex);
      
      if (verbose)
      {
        cout << "Initial cortical parcellation: " << flush;
        cout << seedRegionsTex[0].nItem() << " values, " << flush;
        cout << maxLabel << " regions." << endl;
      }
    }
    else
    {
      throw runtime_error("ERROR: No or wrong seedRegionsTex input.");
    }

    // Read the mask file (if it exists)
    if (roisMaskR.fileName()!= "") roisMaskR.read(roisMask);

    Connectivities * connMatrixToAllMesh_ptr = 0;

    vector<Connectivities* > connMatrixCortexToMesh_ptr_vector;
    connMatrixCortexToMesh_ptr_vector.reserve(meshes_nb);

    if (connMatrixComputingType == "meshintersectionpoint" ||
        connMatrixComputingType == "meshintersectionpointfast")
      connMatrixToAllMesh_ptr = makeConnMatrix_intersection(
          connMatrixComputingType, bundleFilename, inAimsMesh, motionName,
          distthresh, meshClosestPoint_maxDistance, roisMask, length_min,
          length_max, verbose);
    else  // connMatrixComputingType == "meshclosestpoint"
      connMatrixToAllMesh_ptr = makeConnMatrix_closestPoint(
          bundleFilename, inAimsMesh, motionName, verbose);

    if (connectivityTextureType == "seed_connection_density")
      makeConnectivityTexture_seedConnectionDensity(
          connMatrixToAllMesh_ptr, connTextureFileName,
          seedRegionsTex, seedRegionLabel, verbose);
    else if (connectivityTextureType == "seed_mean_connectivity_profile")
      makeConnectivityTexture_seedMeanConnectivityProfile(
          connMatrixToAllMesh_ptr, connTextureFileName, seedRegionsTex,
          seedRegionLabel, maxLabel, distthresh, wthresh, logOption,
          logFile, connMatrixFileName, inAimsMesh, normalize, connMatrixFormat,
          seedRegionVertexIndexType, seedRegionVertexIndexFileName,
          connTextureToTargetMeshesFileNames, meshes_nb,
          connMatrixCortexToMesh_ptr_vector, verbose);

    delete connMatrixToAllMesh_ptr;

    for (int meshLabel = 0; meshLabel < meshes_nb; ++meshLabel)
    {
      // delete connMatrixCortexToMesh_ptr_vector[meshLabel];
    }

    return EXIT_SUCCESS;
  }

  catch (carto::user_interruption &)
  {
    // Exceptions thrown by command line parser (already handled, simply exit)
  }

  catch (exception & e)
  {
    cerr << e.what() << endl;
  }
}

