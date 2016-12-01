#include <constellation/connMatrix.h>
#include <constellation/connMatrixTools.h>
#include <constellation/connConversion.h>
#include <aims/mesh/texturetools.h>
#include <aims/getopt/getopt2.h>
#include <aims/sparsematrix/sparseordensematrix.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

int main(int argc, const char* argv[]) {
  try {
    //-------------------------------------------------------------------------
    Reader< TimeTexture<short> > corticalParcellationReader;
    TimeTexture<short> corticalParcellation;
    int seedLabel = 0;
    int maxLabel  = 0;
    int minLabel  = 0;

    Reader< TimeTexture<short> > secondCorticalParcellationReader;
    TimeTexture<short> secondCorticalParcellation;
    int secondMaxLabel = 0;
    vector<int> secondSeedLabel_vector;
    secondSeedLabel_vector.reserve(1);
    secondSeedLabel_vector.push_back(0);

    string inputMatrixFilename;
    string matrixFormat = "binar_sparse";
    string connectivityType = "seed_connection_density";
    string outputMatrixFilename;
    string outputSeedVertexIndicesTxtFile;
    string outputMeanProfileFilename;
    string outputSeedVertexIndicesTextureFile;

    bool normalize = false;
    bool verbose = false;
    //-------------------------------------------------------------------------
    AimsApplication app(
      argc, argv,
      "Compute the connectivity profile(s). \n\
      Modes supported: (1) mesh_to_mesh (2) oneSeedRegion_to_mesh");

    app.addOption(
      matrixFormat, "-connfmt",
      "Input connectivity matrix format. \n\
      Formats supported: (1) .ima, (2) binar_sparse or (3) ascii_sparse \n\
      default = binar_sparse", true);
    app.addOption(
      inputMatrixFilename, "-connmatrixfile",
      "Input connectivity matrix filename.");
    app.addOption(
      outputMeanProfileFilename, "-outconntex",
      "Output mean connectivity profile filename.");
    app.addOption(
      corticalParcellationReader, "-seedregionstex",
      "Input initial cortical parcellation (e.g. gyri segmentation).");
    app.addOption(
      seedLabel, "-seedlabel",
      "Input seed region label.\n\
      0 to calculate the mean connectivity for all the regions: \n\
      mat.shape = (maxLabel, targetRegionsNb) \n\
      -> if seedLabel != 0 : \n\
      mat.shape = (seedLabel vertexNb, targetRegionsNb)", true);
    app.addOption(
      secondCorticalParcellationReader, "-seedregionstex2",
      "Second cortical parcellation.\n\
      If you want to study a region in particular, \n\
      the input connectivity matrix correspond to the first seed region",
      true);
    app.addOptionSeries(
      secondSeedLabel_vector, "-seedlabel2",
      "second input seed region label + an -seedregiontex2 entry, \n\
      to obtain a connectivity profile of a sub region \n\
      = (of label seedlabel2 in secondCorticalParcellation) \n\
      of the input seedregion (of label seedlabel)", true);
    app.addOption(
      verbose, "-verbose",
      "For more information.", true);
    app.addOption(
      normalize, "-normalize",
      "Normalize connectivity matrix.", true);
    app.addOption(
      connectivityType, "-type",
      "Connectivity types supported: \n\
      (1) seed_mean_connectivity_profile (2) seed_connection_density\n\
      default = seed_connection_density");
    app.addOption(
      outputMatrixFilename, "-outconnmatrix",
      "Output connectivity matrix filename of the seed region,\n\
      size = (seedRegionVertexNb,meshVertexNb)\n\
      default = don\'t save", true);
    app.addOption(
      outputSeedVertexIndicesTxtFile, "-outseedvertex",
      "Output seedlabel region vertex index filename (txt file)", true);
    app.addOption(
      outputSeedVertexIndicesTextureFile, "-outseedvertextex",
      "Output seedlabel region profile filename (texture file)", true);

    app.initialize();
    //-------------------------------------------------------------------------


    // Check the input parameters
    if (matrixFormat != "binar_sparse" && matrixFormat != "ascii_sparse")
      throw runtime_error(string("Unknown matrix format ") + matrixFormat);
    if (connectivityType != "seed_connection_density"
        && connectivityType != "seed_mean_connectivity_profile")
      throw runtime_error(string("Unknown connectivity texture type ")
        + connectivityType);

    // Read the initial cortical parcellation
    corticalParcellationReader.read(corticalParcellation);

    // Define the largest region number
    maxLabel = textureMax(corticalParcellation);

    // Define the smallest region number
    minLabel = textureMin(corticalParcellation);

    // Check that the cortical region is correct
    if (seedLabel < minLabel or seedLabel > maxLabel or seedLabel == 0) {
      throw runtime_error("No or wrong cortical region (seedLabel).");
    }

    // Read the second cortical parcellation, if it exists
    if (secondCorticalParcellationReader.fileName() != "") {
      secondCorticalParcellationReader.read(secondCorticalParcellation);
      secondMaxLabel = textureMax(secondCorticalParcellation);
    }

    // Compute the number of points for each cortical region
    map<short, size_t> *labels_ptr = labelsHistogram(corticalParcellation,
                                                     maxLabel,
                                                     minLabel,
                                                     verbose);
    if (verbose)
      cout << "Number of points of the cortical region: "
        << (*labels_ptr)[seedLabel] << "/" << corticalParcellation[0].nItem()
        << endl;

    // Compute the index of the each cortical region vertex
    Connectivities *connMatrixToAllMesh_ptr = 0;
    vector<size_t> seedVertexIndex;
    seedVertexIndex.reserve((*labels_ptr)[seedLabel]);
    for (size_t i = 0; i < corticalParcellation[0].nItem(); ++i)
      if (corticalParcellation[0][i] == seedLabel)
        seedVertexIndex.push_back(i);
    
    // Read the connectivity matrix
    aims::SparseOrDenseMatrix AllMeshConnMatrix;
    Reader<aims::SparseOrDenseMatrix> spmreader(inputMatrixFilename);
    spmreader.read(AllMeshConnMatrix);

    // If second initial cortical parcellation is specified
    if (secondCorticalParcellationReader.fileName() != "") {
      connMatrixToAllMesh_ptr = new constel::Connectivities(
          (*labels_ptr)[seedLabel],
          til::SparseVector<double>(AllMeshConnMatrix.getSize2()));
      if (AllMeshConnMatrix.getSize1() == (*labels_ptr)[seedLabel]) {
        if (verbose)
          cout << "Connectivity matrix corresponds to seed region" << endl;
        for (size_t i = 0; i < (*labels_ptr)[seedLabel]; ++i) {
          (*connMatrixToAllMesh_ptr)[seedVertexIndex[i]]
            = AllMeshConnMatrix.getSparseRow
              <til::SparseVector<double> >((int32_t)i);
        }
      } else {
        if (verbose)
          cout << "Connectivity matrix corresponds to all the mesh" << endl;
        for (size_t i = 0; i < seedVertexIndex.size(); ++i) {
          (*connMatrixToAllMesh_ptr)[seedVertexIndex[i]]
            = AllMeshConnMatrix.getSparseRow
              <til::SparseVector<double> >((int32_t)seedVertexIndex[i]);
        }
      }
      AllMeshConnMatrix = SparseOrDenseMatrix();
      SparseOrDenseMatrix *mat2 = connectivitiesToSparseOrDenseMatrix(
          *connMatrixToAllMesh_ptr);
      AllMeshConnMatrix = *mat2;
      delete mat2;
    }

    /* Compute the density profile
    - outputDensityTex: texture of the connection density of the seed region
      towards all the other vertex of the mesh
    - outputTargetDensityTex: texture of the connection density of the entire
      mesh towards the seed region */
    if (connectivityType == "seed_connection_density") {
      TimeTexture<float> outputTargetDensityTex;
      outputTargetDensityTex = densityTexture(AllMeshConnMatrix,
                                              seedVertexIndex);
      Writer< TimeTexture<float> > wt(outputMeanProfileFilename);
      wt.write(outputTargetDensityTex);

    } else { // Compute the mean connectivity profile
      // No second cortical parcellation
      if (secondCorticalParcellationReader.fileName().empty()) {
        if (normalize) connMatrixNormalize(AllMeshConnMatrix);
        if (outputMatrixFilename != "") {
          Writer<SparseOrDenseMatrix> w(outputMatrixFilename);
          w.write(AllMeshConnMatrix);
        }

        // Write region vertex indexes (ascii format only)
        if (matrixFormat == "ascii_sparse") {
          if (seedVertexIndex.size() != 0) {
            if (!outputSeedVertexIndicesTxtFile.empty()) {
              fstream findex;
              ostringstream s;
              findex.open(outputSeedVertexIndicesTxtFile.c_str(),
                          fstream::out);
              for (size_t i = 0; i < (*labels_ptr)[seedLabel]; ++i) {
                findex << seedVertexIndex[i] << endl;
              }
            }

            if (!outputSeedVertexIndicesTextureFile.empty()) {
              TimeTexture<unsigned int> seedRegionVertexIndexTex;
              seedRegionVertexIndexTex[0].reserve(seedVertexIndex.size());
              for (size_t vertex = 0; vertex < seedVertexIndex.size();
                   vertex++) {
                unsigned int val = seedVertexIndex[vertex];
                seedRegionVertexIndexTex[0].push_back(seedVertexIndex[vertex]);
              }
              Writer<TimeTexture<unsigned int > > w(
                outputSeedVertexIndicesTextureFile);
              w.write(seedRegionVertexIndexTex);
            }
          }
        }

        TimeTexture<float> outputTargetDensityTex;
        outputTargetDensityTex = meshDensityTexture(AllMeshConnMatrix);
        Writer<TimeTexture<float> > wt(outputMeanProfileFilename);
        wt.write(outputTargetDensityTex);

      } else { // Second cortical parcellation (sub region of the first input seed region)
        int seedLabel2_nb = secondSeedLabel_vector.size();
        TimeTexture<float> outputTargetDensityTex_allLabels2(
          seedLabel2_nb, corticalParcellation[0].nItem());
        for (int i = 0; i < seedLabel2_nb; ++i) {
          int secondSeedLabel = secondSeedLabel_vector[i];

          map<short, size_t> *labels2_ptr = labelsHistogram(
              secondCorticalParcellation,
              secondMaxLabel,
              minLabel,
              verbose);

          vector<size_t> *secondSeedVertexIndex = 0;
          size_t seedRegionLabelVertexNb2 = (*labels2_ptr)[secondSeedLabel];

          Connectivities * extractConnMatrix_ptr
            = connMatrixRegionExtract(
              AllMeshConnMatrix, secondCorticalParcellation, secondSeedLabel,
              seedRegionLabelVertexNb2, & secondSeedVertexIndex);
          AllMeshConnMatrix = SparseOrDenseMatrix();
          SparseOrDenseMatrix *mat2 = connectivitiesToSparseOrDenseMatrix(
              *extractConnMatrix_ptr);
          AllMeshConnMatrix = *mat2;

          delete mat2;
          delete extractConnMatrix_ptr;

          if (normalize) connMatrixNormalize(AllMeshConnMatrix);

          TimeTexture<float> outputTargetDensityTex
            = meshDensityTexture(AllMeshConnMatrix);
          size_t tex_size = outputTargetDensityTex[0].nItem();
          for (size_t v = 0; v < tex_size; ++v) {
            outputTargetDensityTex_allLabels2[i].item(v)
              = outputTargetDensityTex[0][v]/float(seedRegionLabelVertexNb2);
          }
          delete secondSeedVertexIndex;
          delete labels2_ptr;
        }
        Writer<TimeTexture<float> > wt(outputMeanProfileFilename);
        wt.write(outputTargetDensityTex_allLabels2);
      }
    }
    delete labels_ptr;
    delete connMatrixToAllMesh_ptr;

    return EXIT_SUCCESS;
  }
  catch(user_interruption &) {}

  catch(exception & e) {
    cerr << e.what() << endl;
  }

  return EXIT_FAILURE;
}
