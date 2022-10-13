#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <aims/fibers/bundlesFusion.h>
#include <constellation/selectFiberListenerFromMesh.h>
#include <constellation/selectBundlesFromNames.h>
#include <constellation/selectBundlesFromLength.h>
#include <constellation/fibersWeightsReader.h>
#include <constellation/fibersWeightsWriter.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;



int main(int argc, const char* argv[]) {
  try {
    vector<string> fileNameIn;
    string fileNameOut;
    string fileNameOut_notinmesh;
    string weightsFilename;
    string mode = "Name1_Name2orNotInMesh";
    string mname;
    string gyrus;
    rc_ptr<TimeTexture<short> > tex;
    Reader<AimsSurfaceTriangle> r;
    Reader<TimeTexture<short> > texR;
    vector<string> namesList;
    AffineTransformation3d motion;
    int addInt = 0;
    bool as_regex = false;
    bool verbose  = false;
    float cortMinlength =  0;
    float cortMaxlength = -1;
    float nimMinlength  = -1;
    float nimMaxlength  = -2;

    AimsApplication app(
        argc, argv,
        "Selection of tracking fibers according to a white matter mesh with \n\
        label texture.");

    app.addOptionSeries(fileNameIn, "-i", "Bundle File input");
    app.addOption(fileNameOut, "-o", "Bundle File output");
    app.addOption(
        fileNameOut_notinmesh, "-n",
        "Bundle File output for \"distant fibers\"");
    app.addOption(
        r, "--mesh",
        "input mesh");
    app.addOption(
        texR, "--tex",
        "labels input texture" );
    app.addOption(
        weightsFilename, "--weightsFilename",
        "Weights text file matching bundle file input", true);
    app.addOption(
        gyrus, "-g",
        "gyrus name for distant fibers filtering");
    app.addOption(
        mname, "--trs",
        "transformation from t2 to anat", true);
    app.addOption(
        addInt, "--intadd",
        "int to add to Bundle Names", true);
    app.addOption(
        mode, "--mode",
        "mode of new bundles names : = Name1_Name2 or NameFront, or NameEnd,\n\
        or NameFront_NameEnd, or Name1_Name2orNotInMesh (default)", true);
    app.addOptionSeries(
        namesList, "--names",
        "list of string containing bundle names to be selected");
    app.addOption(
        as_regex, "-r",
        "names are regular expressions", true);
    app.addOption(
        cortMinlength, "-l",
        "minimum length for a \"near cortex\" fiber (default: 0)", true);
    app.addOption(
        cortMaxlength, "-L",
        "maximum length for a \"near cortex\" fiber (default: no max)", true);
    app.addOption(
        nimMinlength, "--nimlmin",
        "minimum length for a \"not in mesh\" fiber (default: same as cortex)",
        true);
    app.addOption(
        nimMaxlength, "--nimlmax",
        "maximum length for a \"not in mesh\" fiber (default: same as cortex)",
        true);
    app.addOption(
        verbose, "--verbose",
        "show as much information as possible", true);

    app.initialize();

    rc_ptr<AimsSurfaceTriangle> mesh( r.read() );

    if (nimMinlength < 0) nimMinlength = cortMinlength;
    if (nimMaxlength == -2) nimMaxlength = cortMaxlength;

    if (!mname.empty()) {
      Reader<AffineTransformation3d> mreader(mname);
      mreader.read(motion);
    }

    if (verbose) {
      cout << "reading texture..." << flush;
    }
    tex.reset(texR.read());

    if (verbose) {
      cout << "done" << endl;
      cout << "# vertices : " << mesh->vertex().size() << endl;
      cout << "# faces : " << mesh->polygon().size() << endl;
      cout << "Texture dim : " << (*tex)[0].nItem() << endl;
    }

    /*  Pipeline structure:
        For each input file:

              bundle: BundleReader, reads fibers
                               |
            gyriFilter: SelectFiberListenerFromMesh
                set gyrus pair name on each fiber
                               |
                              / \
                             /   \
                (near cortex)     (not in mesh)
         selectCortexBundles:     selectNimBundles:
       SelectBundlesFromNames     SelectBundlesFromNames
         keeps cortex bundles     keeps bundles going outside
                  |                      |
    selectCBundlesFromLength:     selectNBundlesFromLength:
      SelectBundlesFromLength     SelectBundlesFromLength
             filter by length     filter by length
                  |                      |
 (other files)    |                      |   (other files)
         \        |                      |    /
            cortexRegroup:        nimRegroup:
             BundleFusion         BundleFusion
         sorts fibers into        sorts fibers into
        consistent bundles        consistent bundles
                  |                      |
            cortexWriter          nimWriter
    */

    unsigned i, n = fileNameIn.size();

    // MRtrix pipeline
    // Read fiber weight if file is provided
    if ( !weightsFilename.empty() &&
         fileNameIn[0].substr(fileNameIn[0].size() - 4) == ".tck") {
        // TODO: handle multiple files given

        // create common elements
        rc_ptr<BundleReader> bundle;
        rc_ptr<SelectFiberListenerFromMesh> gyriFilter;
        rc_ptr<FibersWeightsReader> fibersWeightsReader;

        // Cortex
        rc_ptr<SelectBundlesFromNames> selectCortexBundles;
        rc_ptr<SelectBundlesFromLength> selectCBundlesFromLength;
        rc_ptr<BundlesFusion> cortexRegroup(new BundlesFusion(1));
        rc_ptr< FibersWeightsWriter > cortexWeightsWriter;
        rc_ptr< BundleWriter > cortexWriter(new BundleWriter);

        // NIM
        rc_ptr<SelectBundlesFromNames> selectNimBundles;
        rc_ptr<SelectBundlesFromLength> selectNBundlesFromLength;
        rc_ptr<BundlesFusion> nimRegroup(new BundlesFusion(1));
        rc_ptr< FibersWeightsWriter > nimWeightsWriter;
        rc_ptr< BundleWriter > nimWriter(new BundleWriter);

        // rc_ptrs must be allocated outside of any block, otherwise they will
        // get deleted at the end if the block (if() { }):
        // bundle producers only keep pointers, not rc_ptrs thus do not maintain
        // life of their listeners

        // Bundles reader creation
        string fileName = fileNameIn[0];
        bundle.reset( new BundleReader( fileName ) );

        // Weights reader
        fibersWeightsReader.reset( new FibersWeightsReader( weightsFilename ));
        bundle->addBundleListener(*fibersWeightsReader);

        // Set names from mesh/label texture
        gyriFilter.reset(
          new SelectFiberListenerFromMesh(mesh, tex, mode, addInt, motion,""));
        fibersWeightsReader->addBundleListener(*gyriFilter);

        // -- 1st branch: near cortex
        // filter labels
        selectCortexBundles.reset(
          new SelectBundlesFromNames(namesList, verbose, as_regex, true));
        gyriFilter->addBundleListener(*selectCortexBundles);

        // filter from length
        selectCBundlesFromLength.reset(
          new SelectBundlesFromLength(cortMinlength, cortMaxlength, verbose));
        selectCortexBundles->addBundleListener(*selectCBundlesFromLength);

        // cortex bundle fusion
        selectCBundlesFromLength->addBundleListener(*cortexRegroup);

        // cortex fibers weights writer
        string cortexWeightsFilename = fileNameOut.substr(
          0, fileNameOut.size() - 8) + "_weights.txt";
        cortexWeightsWriter.reset(
          new FibersWeightsWriter(cortexWeightsFilename));
        cortexRegroup->addBundleListener(*cortexWeightsWriter);

        // cortex bundles writer
        cortexWriter->setFileString(fileNameOut);
        cortexWeightsWriter->addBundleListener(*cortexWriter);

        // -- 2nd branch: not in mesh
        // filter labels
        vector<string> notinmesh_names;
        notinmesh_names.push_back(gyrus + "_notInMesh");
        selectNimBundles.reset(
          new SelectBundlesFromNames(notinmesh_names, verbose, false, true));
        gyriFilter->addBundleListener(*selectNimBundles);

        // filter from length
        selectNBundlesFromLength.reset(
          new SelectBundlesFromLength(nimMinlength, nimMaxlength, verbose));
        selectNimBundles->addBundleListener(*selectNBundlesFromLength);

        // nim bundle fusion
        selectNBundlesFromLength->addBundleListener(*nimRegroup);

        // NIM fibers weights writer
        string nimWeightsFilename = fileNameOut_notinmesh.substr(
          0, fileNameOut_notinmesh.size() - 8) + "_weights.txt";
        nimWeightsWriter.reset(new FibersWeightsWriter(nimWeightsFilename));
        nimRegroup->addBundleListener(*nimWeightsWriter);

        // NIM bundles writer
        nimWriter->setFileString(fileNameOut_notinmesh);
        nimWeightsWriter->addBundleListener(*nimWriter);

        // Run pipeline
        if (verbose) cout << "process file: " << fileName << endl;
        bundle->read();
        if (verbose) cout << "done for " << fileName << endl;
        if (verbose) cout << "MRtrix pipeline done.\n";
    }
    // Connectomist pipeline
    else
    {
      // create common elements
      vector<rc_ptr<BundleReader> > bundle(n);
      vector<rc_ptr<SelectFiberListenerFromMesh> > gyriFilter(n);

      // cortex
      vector<rc_ptr<SelectBundlesFromNames> > selectCortexBundles(n);
      vector<rc_ptr<SelectBundlesFromLength> > selectCBundlesFromLength(n);
      // regroup cortex bundles by gyrus name
      rc_ptr<BundlesFusion> cortexRegroup(new BundlesFusion((int) n));
      rc_ptr< BundleWriter > cortexWriter(new BundleWriter);

      // nim
      vector<rc_ptr<SelectBundlesFromNames> > selectNimBundles(n);
      vector<rc_ptr<SelectBundlesFromLength> > selectNBundlesFromLength(n);
      // regroup "not in mesh" bundles by gyrus name
      rc_ptr<BundlesFusion> nimRegroup( new BundlesFusion((int) n));
      rc_ptr< BundleWriter > nimWriter(new BundleWriter);

      // duplicate branches for all input files
      for (i=0; i!=n; ++i) {
        // Bundles reader creation
        string fileName = fileNameIn[i];
        bundle[i].reset( new BundleReader( fileName ) );

        // Set names from mesh/label texture
        gyriFilter[i].reset(
          new SelectFiberListenerFromMesh(mesh, tex, mode, addInt, motion,""));
        bundle[i]->addBundleListener(*gyriFilter[i]);

        // -- 1st branch: near cortex
        // filter labels
        selectCortexBundles[i].reset(
          new SelectBundlesFromNames(namesList, verbose, as_regex, true));
        gyriFilter[i]->addBundleListener(*selectCortexBundles[i]);

        // filter from length
        selectCBundlesFromLength[i].reset(
          new SelectBundlesFromLength(cortMinlength, cortMaxlength, verbose));
        selectCortexBundles[i]->addBundleListener(*selectCBundlesFromLength[i]);

        // connect to common cortex fusion element
        selectCBundlesFromLength[i]->addBundleListener(*cortexRegroup);

        // -- 2nd branch: not in mesh
        // filter labels
        vector<string> notinmesh_names;
        notinmesh_names.push_back(gyrus + "_notInMesh");
        selectNimBundles[i].reset(
          new SelectBundlesFromNames(notinmesh_names, verbose, false, true));
        gyriFilter[i]->addBundleListener(*selectNimBundles[i]);

        // filter from length
        selectNBundlesFromLength[i].reset(
          new SelectBundlesFromLength(nimMinlength, nimMaxlength, verbose));
        selectNimBundles[i]->addBundleListener(*selectNBundlesFromLength[i]);

        // regroup bundles by gyrus name
        selectNBundlesFromLength[i]->addBundleListener(*nimRegroup);
      }
      // cortex bundles writer
      cortexWriter->setFileString(fileNameOut);
      cortexRegroup->addBundleListener(*cortexWriter);

      // NIM bundles writer
      nimWriter->setFileString(fileNameOut_notinmesh);
      nimRegroup->addBundleListener(*nimWriter);

      // run all those
      for (i = 0; i < n; ++i) {
        if (verbose) cout << "process file: " << fileNameIn[i] << endl;
        bundle[i]->read();
        if (verbose) cout << "done for " << fileNameIn[i] << endl;
      }

      if (verbose) cout << "Connectomist pipeline done.\n";
    }

  return EXIT_SUCCESS;
  }
  catch (carto::user_interruption &) {
    // Exceptions thrown by command line parser (already handled, simply exit)
  }
  catch (exception & e) {
    cerr << e.what() << endl;
  }
}
