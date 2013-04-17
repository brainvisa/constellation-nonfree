#include <constellation/connMatrix.h>
#include <constellation/connMatrixTools.h>
#include <aims/mesh/texturetools.h>
#include <aims/getopt/getopt2.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

int main( int argc, const char* argv[] )
{
  try
  {
    Reader<AimsSurfaceTriangle> inMeshAimsR;
    Reader< TimeTexture<short> > seedRegionsTexR;
    AimsSurfaceTriangle inAimsMesh;
    TimeTexture<short> seedRegionsTex;
    size_t seedRegionLabel = 0;

    size_t seedRegionsNb = 0;
    string connTextureFileName;
    string outConnMatrixFileName = "";
    bool verbose = false;
    bool normalize = false;
    string connectivityTextureType = "seed_connection_density";
    string outSeedRegionVertexIndexFileName = "";
    string outSeedRegionVertexIndexFileNameTex = "";
    string connMatrixFormat = "binar_sparse";
    string inConnMatrixFileName;

    Reader< TimeTexture<short> > seedRegionsTexR2;
    TimeTexture<short> seedRegionsTex2;
    vector<int> seedRegionLabel2_vector;
    seedRegionLabel2_vector.reserve(1);
    seedRegionLabel2_vector.push_back(0);
    size_t seedRegionsNb2 = 0;

    float distthresh = 0;
    string outSmoothedConnMatrixFileName = "";

    AimsApplication app( argc, argv,
      "This command computes and writes connection density texture(s). There are two modes: mesh_to_mesh, oneSeedRegion_to_mesh" );
    app.addOption( inMeshAimsR, "-mesh", "input mesh" );
    app.addOption( connMatrixFormat, "-connfmt",
      "input conn matrix format, .ima or binar_sparse or ascii_sparse; default = binar_sparse",true );
    app.addOption( inConnMatrixFileName, "-connmatrixfile",
      "input conn matrix filename" );
    app.addOption( connTextureFileName, "-outconntex",
      "output mean connectivity texture file name" );
    app.addOption( seedRegionsTexR, "-seedregionstex",
      "input region texture" );
    app.addOption( seedRegionLabel, "-seedlabel", "input seed region label, 0 to calculate the mean connectivity for all the regions : connMatrix.shape = (seedRegionsNb, targetRegionsNb); = (seedRegionLabel vertex nb, targetRegionsNb) if seed region label != 0", true );
    app.addOption( seedRegionsTexR2, "-seedregionstex2", "second input region texture : if you want to study a region in particular, in case of the input conn matrix correspond to the first seed region (conn matrix shape = (seedRegion Of SeedLabel seedlabel vertex nb, mesh vertex nb", true );
    app.addOptionSeries( seedRegionLabel2_vector, "-seedlabel2", "second input seed region label, in case of an -seedregiontex2 entry, to obtain a connectivity profile of a sub region (of lable seedlabel2 in seedregionstex2) of the input seedregion (of label seedlabel)", true );
    app.addOption( verbose, "-verbose", "show as much information as possible",
      true );
    app.addOption( normalize, "-normalize", "normalize connectivity matrix",
      true);
    app.addOption( connectivityTextureType, "-type",
      "connectivity type: seed_mean_connectivity_profile or seed_connection_density default = seed_connection_density" );
    app.addOption( outConnMatrixFileName, "-outconnmatrix",
      "output connectivity matrix filename of the seed region, size = (seedRegionVertexNb,meshVertexNb), default = don t save", true );
    app.addOption( outSeedRegionVertexIndexFileName, "-outseedvertex",
      "output seedlabel region vertex index file name (txt file)", true );
    app.addOption( outSeedRegionVertexIndexFileNameTex, "-outseedvertextex",
      "output seedlabel region texture file name (texture file)", true );
    app.initialize();

    if( connMatrixFormat != "binar_sparse"
        && connMatrixFormat != "ascii_sparse" )
      throw runtime_error( string( "unknown matrix format " )
        + connMatrixFormat );
    if( connectivityTextureType != "seed_connection_density"
        && connectivityTextureType != "seed_mean_connectivity_profile" )
      throw runtime_error( string( "unknown connectivity texture type " )
        + connectivityTextureType );

    //Reading inputs
    if(verbose) cout << "Reading inputs :" << endl;

    inMeshAimsR.read( inAimsMesh );
    if(verbose)
      cout << "cortex mesh vertex nb: " << int(inAimsMesh.vertex().size())
        << endl;

    uint cortexMeshVertexNb = uint(inAimsMesh.vertex().size());

    if(verbose)
      cout << "reading seedRegionsTex: " << flush;
    seedRegionsTexR.read( seedRegionsTex );
    if(verbose)
      cout << "Texture dim : " << seedRegionsTex[0].nItem()
        << flush;
    seedRegionsNb = textureMax(seedRegionsTex, inAimsMesh);
    if(verbose)
    {
      cout << ", seedRegionsNb:" << seedRegionsNb << flush;
      cout << "...done." << endl;
    }

    if(seedRegionsTexR2.fileName()!="")
    {
      if(verbose)
        cout << "reading seedRegionsTex2: " << flush;
      seedRegionsTexR2.read( seedRegionsTex2 );
      if(verbose)
        cout << "Texture dim : " << seedRegionsTex2[0].nItem() << flush;
      seedRegionsNb2 = textureMax(seedRegionsTex2, inAimsMesh);
      if(verbose)
        cout << ", seedRegionsNb2:" << seedRegionsNb2 << flush;
      if(verbose)
        cout << "...done." << endl;
    }
    else
    {
      if(verbose)
        cout << "no seedRegionsTex2 input...!" << endl;
    }

    if (verbose)
    {
      cout << "seedRegionLabel : " << seedRegionLabel << endl;
      cout << "Reading inputs done." << endl << endl;
    }

    vector< size_t > * labels_ptr
      = labelsHistogram( seedRegionsTex, seedRegionsNb, verbose );
    vector< size_t > & labels = *labels_ptr;
    Connectivities * connMatrixToAllMesh_ptr = 0;

    vector<size_t> seedVertexIndex;
    seedVertexIndex.reserve( labels[seedRegionLabel] );
    if(verbose)
      cout << "labels[seedRegionLabel]:" << labels[seedRegionLabel] << endl;
    for(size_t i = 0; i < seedRegionsTex[0].nItem(); ++i)
    {
      if (seedRegionsTex[0][i]==seedRegionLabel)
      {
        seedVertexIndex.push_back(i);
      }
    }

    if( verbose )
      cout << "reading matrix...\n";
    aims::SparseMatrix AllMeshConnMatrix;
    if( connMatrixFormat == "binar_sparse" )
      AllMeshConnMatrix.read( inConnMatrixFileName );
    else
      AllMeshConnMatrix.read( inConnMatrixFileName, "ascii" );

    connMatrixToAllMesh_ptr
      = new constel::Connectivities( labels[seedRegionLabel],
          til::SparseVector<double>(AllMeshConnMatrix.getSize2() ) );

    if( seedRegionsTexR2.fileName().empty() )
    {
      Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
      if( AllMeshConnMatrix.getSize1() == labels[seedRegionLabel] )
      {
        if (verbose)
          cout << "input connectivity matrix corresponds to seed region"
            << endl;
        for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
        {
          connMatrixToAllMesh[i] = AllMeshConnMatrix.getSparseRow((int32_t)i);
        }
      }
      else
      {
        if (verbose)
          cout << "input connectivity matrix corresponds to all the mesh"
          << endl;
        for (size_t i = 0; i < seedVertexIndex.size(); ++i)
        {
          connMatrixToAllMesh[i]
            = AllMeshConnMatrix.getSparseRow((int32_t)seedVertexIndex[i]);
        }
      }
    }
    else // seedRegionsTexR2 is specified
    {
      Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
      if(AllMeshConnMatrix.getSize1() == labels[seedRegionLabel])
      {
        if (verbose)
          cout << "input connectivity matrix corresponds to seed region"
            << endl;
        for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
        {
          connMatrixToAllMesh[seedVertexIndex[i]]
            = AllMeshConnMatrix.getSparseRow((int32_t)i);
        }
      }
      else
      {
        if (verbose)
          cout << "input connectivity matrix corresponds to all the mesh"
            << endl;
        for (size_t i = 0; i < seedVertexIndex.size(); ++i)
        {
          connMatrixToAllMesh[seedVertexIndex[i]]
            = AllMeshConnMatrix.getSparseRow((int32_t)seedVertexIndex[i]);
        }
      }
    }
    if(verbose)
      cout << "OK" << endl;


    if (connectivityTextureType == "seed_connection_density")
    {
      /*
      Computing densityTexture:
      */
      TimeTexture<float> outputTargetDensityTex
        = densityTexture( connMatrixToAllMesh_ptr, seedVertexIndex );
      if (verbose)
        cout << "Writing connectivity Density texture:" << connTextureFileName
          << endl;
      Writer< TimeTexture<float> > wt( connTextureFileName );
      wt.write( outputTargetDensityTex );
    }
    else // connectivityTextureType == "seed_mean_connectivity_profile"
    {

      if (seedRegionLabel <= 0 or seedRegionLabel > seedRegionsNb)
      {
        throw runtime_error("no or wrong seedRegionLabel...!");
      }
      //Computing mat [region vertex][meshVertexNb]
      // Read a regions (gyri for example) labeled texture and make a Label Histogram
      // Assumes the regions (gyri) image is labeled between 1 and regionsNb, background = 0 or lower than 1
      // Count all the labels greater than regionsNb and the labels of the background

      vector< size_t > * labels_ptr
        = labelsHistogram( seedRegionsTex, seedRegionsNb, verbose);
      vector< size_t > & labels = *labels_ptr;

      if( seedRegionsTexR2.fileName().empty() ) //no second input region
      {
        size_t seedRegionLabelVertexNb = labels[seedRegionLabel];
        if (normalize)
        {
          Connectivities * connN 
            = connMatrixNormalize(connMatrixToAllMesh_ptr);
          delete connMatrixToAllMesh_ptr;
          connMatrixToAllMesh_ptr = connN;
        }
        if( outConnMatrixFileName != "" )
        {
          if(verbose)
            cout << "Writing seed region connMatrix as binar sparse matrix format:" << outConnMatrixFileName << endl;
          int32_t size1 = (*connMatrixToAllMesh_ptr).size();
          int32_t size2 = (*connMatrixToAllMesh_ptr)[0].size();
          Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
          if(verbose)
            cout << "size1: " << size1 << ", size2: " << size2 << endl;
          SparseMatrix *m
            = connectivitiesToSparseMatrix( *connMatrixToAllMesh_ptr );
          m->write(outConnMatrixFileName);
          delete m;
        }

        // Write region vertex indexes (ascii format only)
        size_t seedVertexIndex_size = seedVertexIndex.size();
        if (seedVertexIndex.size()!=0)
        {
          if( !outSeedRegionVertexIndexFileName.empty() )
          {
            fstream findex;
            ostringstream s;
            findex.open( outSeedRegionVertexIndexFileName.c_str(),
                         fstream::out );
            for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
            {
              findex << seedVertexIndex[i] << endl;
            }
          }
          if( !outSeedRegionVertexIndexFileNameTex.empty() )
          {
            TimeTexture< unsigned int > seedRegionVertexIndexTex;
            cout << "seedVertexIndex_size:" << seedVertexIndex_size << endl;
            seedRegionVertexIndexTex[0].reserve(seedVertexIndex_size);
            for (size_t vertex = 0; vertex < seedVertexIndex_size; vertex++)
            {
              unsigned int val = seedVertexIndex[vertex];
              seedRegionVertexIndexTex[0].push_back(seedVertexIndex[vertex]);
            }
            cout << "seedRegionVertexIndexTex0.nItem():"
              << seedRegionVertexIndexTex.nItem() << ", "
              << seedRegionVertexIndexTex[0].nItem() << endl;
            Writer<TimeTexture<unsigned int > > w(
              outSeedRegionVertexIndexFileNameTex );
            w.write(seedRegionVertexIndexTex);
          }
        }
        else
        {
          cout << "seedVertexIndex size = 0" << endl;
        }

        /*
        Computing densityTexture:
        outputDensityTex: texture of the connection density of the seed region towards all the other vertex of the mesh
        outputTargetDensityTex: texture of the connection density of the entire mesh towards the seed region
        */
        TimeTexture<float> outputTargetDensityTex
          = meshDensityTexture( connMatrixToAllMesh_ptr );
        if (verbose)
          cout << "Writing mean connectivity profile texture:"
            << connTextureFileName << endl;
        Writer<TimeTexture<float> > wt( connTextureFileName );
        wt.write( outputTargetDensityTex );

      }
      else //there is a second input region (sub region of the first input seed region)
      {

        int seedLabel2_nb = seedRegionLabel2_vector.size();
        TimeTexture<float> outputTargetDensityTex_allLabels2(
          seedLabel2_nb, cortexMeshVertexNb);
        for (int i = 0; i < seedLabel2_nb; ++i)
        {
          int seedRegionLabel2 = seedRegionLabel2_vector[i];

          if(verbose)
          {
            cout << "seedRegionLabel2:" << seedRegionLabel2 << endl;
            cout << "seedRegionsNb2:" << seedRegionsNb2 << endl;
          }
          vector< size_t > * labels2_ptr
            = labelsHistogram(seedRegionsTex2, seedRegionsNb2, verbose);
          vector< size_t > & labels2 = *labels2_ptr;

          vector<size_t> * seedVertexIndex2 = 0;
          size_t seedRegionLabelVertexNb2 = labels2[seedRegionLabel2];
          if(verbose)
            cout << "seedRegionLabelVertexNb2:" << seedRegionLabelVertexNb2
              << endl;
          Connectivities * extractConnMatrix_ptr
            = connMatrixRegionExtract(
              connMatrixToAllMesh_ptr, seedRegionsTex2, seedRegionLabel2,
              seedRegionLabelVertexNb2, & seedVertexIndex2 );
          if(normalize)
            extractConnMatrix_ptr = connMatrixNormalize(extractConnMatrix_ptr);
          TimeTexture<float> outputTargetDensityTex
            = meshDensityTexture(extractConnMatrix_ptr);
          if (verbose)
            cout << "Writing mean connectivity profile texture:"
              << connTextureFileName << endl;
          size_t tex_size = outputTargetDensityTex[0].nItem();
          cout << "tex_size:" << tex_size << endl;
          for (size_t v = 0; v < tex_size; ++v)
          {
            outputTargetDensityTex_allLabels2[i].item(v)
              = outputTargetDensityTex[0][v]/float(seedRegionLabelVertexNb2);
          }
          delete seedVertexIndex2;
          delete labels2_ptr;
          delete extractConnMatrix_ptr;
        }
        Writer<TimeTexture<float> > wt( connTextureFileName );
        wt.write( outputTargetDensityTex_allLabels2 );
      }
      delete labels_ptr;
    }

    delete connMatrixToAllMesh_ptr;

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

