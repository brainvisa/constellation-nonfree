#include <constellation/connMatrix.h>
#include <aims/sparsematrix/sparseMatrix.h>
#include <aims/mesh/texturetools.h>
#include <aims/getopt/getopt2.h>

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

    Reader<AimsSurfaceTriangle> inMeshAimsR;
    Reader< TimeTexture<short> > seedRegionsTexR;
    Reader< TimeTexture<short> > targetRegionsTexR;
    AimsSurfaceTriangle inAimsMesh;
    TimeTexture<short> seedRegionsTex;
    TimeTexture<short> targetRegionsTex;
    std::size_t targetRegionsNb = 0;
    std::size_t seedRegionLabel = 0;
    std::size_t seedRegionsNb = 0;
    std::string connTextureFileName;
    std::string connMatrixFileName = "";
    bool verbose = false;
    bool normalize = false;
    std::string connectivityTextureType = "oneSeedRegion_to_targets";
    std::string connLineTargetsDensityFileName = "";
    std::string seedRegionVertexIndexFileName = "";
    std::string connMatrixFormat = "binar_sparse";
    std::string inConnMatrixFileName = "";

    AimsApplication app( argc, argv, "This command computes and writes connection density texture(s). There are two modes: oneSeedRegion_to_targets, seedVertex_to_targets" );
    app.addOption( inMeshAimsR, "-mesh", "input mesh" );
    app.addOption( inConnMatrixFileName, "-connmatrixfile", "input conn matrix filename" );
    app.addOption( seedRegionsTexR, "-seedregionstex", "input region texture : study a region in particular" );
    app.addOption( connTextureFileName, "-outconntex", "output mean connectivity texture file name", true );
    app.addOption( connMatrixFormat, "-connfmt", "input conn matrix format, .ima, binar_sparse, or ascii_sparse; default = binar_sparse", true );
    app.addOption( seedRegionLabel, "-seedlabel", "input seed region label, 0 to calculate the mean connectivity for all the regions : connMatrix.shape = (seedRegionsNb, targetRegionsNb); = (seedRegionLabel vertex nb, targetRegionsNb) if seed region label != 0", true );
    app.addOption( verbose, "-verbose", "show as much information as possible", true );
    app.addOption( normalize, "-normalize", "normalize connectivity matrix", true);
    app.addOption( connectivityTextureType, "-type", "connectivity type: oneSeedRegion_to_targets, seedVertex_to_targets, default = oneSeedRegion_to_targets", true);
    app.addOption( targetRegionsTexR, "-targetregionstex", "target regions texture for the calculation of the global connectivity matrix", true );
    app.addOption( connLineTargetsDensityFileName, "-outconntargets", "output connectivity targets texture file name (shape (1, targetsNb))", true );
    app.addOption( seedRegionVertexIndexFileName, "-seedvertexindex", "seedlabel region vertex indexes file name (text format)", true );
    app.addOption( connMatrixFileName, "-connmatrix", "output connectivity Matrix file name", true );
    app.initialize();

    if( connMatrixFormat != ".ima" && connMatrixFormat != "binar_sparse"
      && connMatrixFormat != "ascii_sparse" )
      throw runtime_error( "incorrect connfmt parameter" );

    //Reading inputs
    if(verbose) std::cout << "Reading inputs :" << std::endl;

    inMeshAimsR.read( inAimsMesh );
    if(verbose) std::cout << "cortex mesh vertex nb:" << int(inAimsMesh.vertex().size()) << std::endl;
    
    if(verbose) std::cout << "reading seedRegionsTex: " << flush;
    seedRegionsTexR.read( seedRegionsTex );
    if(verbose) std::cout << "Texture dim : " << seedRegionsTex[0].nItem() << std::flush;
    seedRegionsNb = textureMax(seedRegionsTex);
    if(verbose) std::cout << ", seedRegionsNb:" << seedRegionsNb << std::flush;
    if(verbose) std::cout << "...done." << std::endl;
    
    std::vector< std::size_t > * labels_ptr = labelsHistogram(seedRegionsTex, seedRegionsNb, verbose);
    std::vector< std::size_t > & labels = *labels_ptr;
    
    Matrix inConnMatrixImage;
    constel::Connectivities * connMatrixToAllMesh_ptr = 0;
    
    std::vector<size_t> seedVertexIndex;
    seedVertexIndex.reserve(labels[seedRegionLabel]);
    if(verbose) std::cout << "labels[seedRegionLabel]:" << labels[seedRegionLabel] << std::endl;
    
    
    if (seedRegionVertexIndexFileName !="")
    {
      std::ifstream _file(seedRegionVertexIndexFileName.c_str());
      std::string vertex_index;
      uint vertexIndex_count = 0;
      while(_file >> vertex_index)
      {
        uint vertexIndex_uint = strtoul(vertex_index.c_str(), NULL, 0);
        seedVertexIndex.push_back(vertexIndex_uint);
//         if(verbose) std::cout << vertex_index << std::endl;
      }
    }
    else
    {
      for(size_t i = 0; i < seedRegionsTex[0].nItem(); ++i)
      {
        if (seedRegionsTex[0][i]==seedRegionLabel)
        {
          seedVertexIndex.push_back(i);
        }
      }
    }
    std::cout << "seedVertexIndex.size():" << seedVertexIndex.size() << std::endl;
    if (seedVertexIndex.size()!=labels[seedRegionLabel])
    {
      throw runtime_error("seedVertexIndex.size()!=labels[seedRegionLabel]");
    }

    uint cortexMeshVertexNb = uint(inAimsMesh.vertex().size());

    if( connMatrixFormat == ".ima" )
    {
      Reader< Matrix > inConnMatrixR( inConnMatrixFileName );
      std::cout << "connMatrixFormat " << connMatrixFormat << std::endl;
      inConnMatrixR.read( inConnMatrixImage );
      if(verbose) std::cout << "conn matrix im dim:" << inConnMatrixImage.dimX() <<", " << inConnMatrixImage.dimY() << ", " << inConnMatrixImage.dimZ() << std::endl;
      if (inConnMatrixImage.dimX() != labels[seedRegionLabel] or inConnMatrixImage.dimY() != cortexMeshVertexNb)
      {
        throw runtime_error("inConnMatrixImage.dimX() != labels[seedRegionLabel] or inConnMatrixImage.dimY() != cortexMeshVertexNb");
      }
      connMatrixToAllMesh_ptr = new constel::Connectivities(inConnMatrixImage.dimX(), til::SparseVector<double>(inConnMatrixImage.dimY()));
      Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
      for(int x=0;x<inConnMatrixImage.dimX();++x)
      {
        for(int y=0;y<inConnMatrixImage.dimY();++y)
        {
          double val = double(inConnMatrixImage(x,y,0,0));
          if (val!=0)
          {
            connMatrixToAllMesh[x][y] = val;
          }
        }
      }
    }
    else
    {
      // sparse matrix formats
      std::cout << "sparse format for connectivity matrix:"
        << inConnMatrixFileName << std::endl;
      std::cout << "seedRegionVertexIndexFileName:"
        << seedRegionVertexIndexFileName << std::endl;
      aims::SparseMatrix AllMeshConnMatrix;
      if (connMatrixFormat == "binar_sparse")
      {
        AllMeshConnMatrix.read(inConnMatrixFileName);
      }
      else if (connMatrixFormat == "ascii_sparse")
      {
        AllMeshConnMatrix.read(inConnMatrixFileName, "ascii");
      }

      if(verbose) std::cout << "size1: " << size_t(AllMeshConnMatrix.getSize1()) <<  ", size2: " << size_t(AllMeshConnMatrix.getSize2()) << std::endl;
      connMatrixToAllMesh_ptr = new constel::Connectivities(labels[seedRegionLabel], til::SparseVector<double>(AllMeshConnMatrix.getSize2()));
      Connectivities & connMatrixToAllMesh = *connMatrixToAllMesh_ptr;
      if(AllMeshConnMatrix.getSize1() == labels[seedRegionLabel])
      {
        if (verbose)  std::cout << "input connectivity matrix correspond to seed region" << std::endl;
        for (size_t i = 0; i < labels[seedRegionLabel]; ++i)
        {
          connMatrixToAllMesh[i] = AllMeshConnMatrix.getSparseRow
            <til::SparseVector<double> >((int32_t)i);
        }

      }
      else
      {
        if (verbose)  std::cout << "input connectivity matrix correspond to all the mesh" << std::endl;
        for (size_t i = 0; i < seedVertexIndex.size(); ++i)
        {
//           std::cout << "i:" << i << std::endl;
          connMatrixToAllMesh[i] = AllMeshConnMatrix.getSparseRow
            <til::SparseVector<double> >((int32_t)seedVertexIndex[i]);
        }
      }
    }

    if(targetRegionsTexR.fileName()!="")
    {
      if(verbose) std::cout << "reading targetRegionsTex..." << std::flush;
      targetRegionsTexR.read( targetRegionsTex );
      if(verbose) std::cout << "Texture dim : " << targetRegionsTex[0].nItem() << std::flush;
      targetRegionsNb = textureMax(targetRegionsTex);
      if(verbose) std::cout << ", targetRegionsNb: " << targetRegionsNb << std::flush;
      if(verbose) std::cout << "...done." << std::endl;
    }
    else
    {
      throw runtime_error( "no or wrong targetRegionsTex input...!" );
    }

    if (verbose)
    {
      std::cout << "seedRegionLabel : " << seedRegionLabel << std::endl;
    }
    if(verbose) std::cout << "Reading inputs done." << std::endl << std::endl;
    
    
    if (connectivityTextureType == "oneSeedRegion_to_targets" or connectivityTextureType == "seedVertex_to_targets")
    {
      /*
      Computing densityTexture:
      */
      if(connMatrixToAllMesh_ptr)
      {
        std::cout << "exist:" << std::endl;
      }
      else
      {
        std::cout << "not exist !" << std::endl;
      }
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
      Connectivities * extractAndRegroupConnMatrix_ptr = connMatrixSeedMesh_to_targetMeshTargets_regroup(connMatrixToAllMesh_ptr, targetRegionsTex, targetRegionsNb);
      
      if (connectivityTextureType == "oneSeedRegion_to_targets")
      {
        constel::Connectivity * totalConnSeedRegionToTargets_ptr = connMatrixSumRows(extractAndRegroupConnMatrix_ptr);
        
        if (normalize)
        {
          constel::Connectivity & line_i = *totalConnSeedRegionToTargets_ptr;
          if (line_i.is_null()!=1)
          {
            constel::Connectivity sparseVectorIdentity = constel::Connectivity(line_i.size());
            for(std::size_t j = 0; j < line_i.size(); ++j)
            {
              sparseVectorIdentity[j]=1.;
            }
            double sum_i = til::dot<double>(sparseVectorIdentity, line_i);
            if(sum_i!=0)
            {
              for(std::size_t j = 0; j < line_i.size(); ++j)
              {
                double connval = line_i[j];
                line_i[j] = (1/sum_i)*connval;
              }
            }
          }
        }
        constel::Connectivity & totalConnSeedRegionToTargets = *totalConnSeedRegionToTargets_ptr;
        TimeTexture<float> * outputTargetDensityTex_ptr = oneTargetDensityTargetsRegroupTexture(totalConnSeedRegionToTargets_ptr, targetRegionsTex);
        TimeTexture<float> & outputTargetDensityTex = *outputTargetDensityTex_ptr;
        if (verbose) std::cout << "Writing Density texture, shape (1, meshVertexNb):" << connTextureFileName << std::endl;
        Writer<TimeTexture<float> > wt( connTextureFileName );
        wt.write( outputTargetDensityTex );
        if (connLineTargetsDensityFileName!="")
        {
          TimeTexture<float> targetsDensity_tex;
          targetsDensity_tex[0].reserve(targetRegionsNb);
          for (std::size_t t = 0; t < targetRegionsNb; ++t)
          {
            targetsDensity_tex[0].push_back(totalConnSeedRegionToTargets[t]);
          }
          Writer<TimeTexture<float> > wt2( connLineTargetsDensityFileName );
          wt2.write( targetsDensity_tex );
          if(verbose) std::cout << "targets density connectivity texture written...! " << std::endl;
        }
        else
        {
          if(verbose) std::cout << "Don't write connectivity of targets texture, shape (1, targetsNb)." << std::endl;
        }
        if (connMatrixFileName!= "")
        {
          writeAimsFmtConnMatrix(extractAndRegroupConnMatrix_ptr,connMatrixFileName);
        }
        else
        {
          if(verbose) std::cout << "Don't write connectivity matrix image." << std::endl;
        }
        delete totalConnSeedRegionToTargets_ptr;
        delete outputTargetDensityTex_ptr;
        
      }
      else // (connectivityTextureType == "seedVertex_to_targets")
      {
        if (normalize) extractAndRegroupConnMatrix_ptr = connMatrixNormalize(extractAndRegroupConnMatrix_ptr);
        std::cout << "seedVertex_to_targets connMatrix computed." << std::endl;
        if (connMatrixFileName!= "")
        {
          writeAimsFmtConnMatrix(extractAndRegroupConnMatrix_ptr,connMatrixFileName);
        }
        else
        {
          if(verbose) std::cout << "Don't write connectivity matrix image." << std::endl;
        }
      }
      delete extractAndRegroupConnMatrix_ptr;
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

