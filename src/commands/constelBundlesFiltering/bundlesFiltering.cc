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
      tex.reset(texR.read());
    }

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

    // create common elements

    vector<rc_ptr<BundleReader> > bundle(n);
    vector<rc_ptr<SelectFiberListenerFromMesh> > gyriFilter(n);

    vector<rc_ptr<SelectBundlesFromNames> > selectCortexBundles(n);
    vector<rc_ptr<SelectBundlesFromLength> > selectCBundlesFromLength(n);
    // regroup cortex bundles by gyrus name
    rc_ptr<BundlesFusion> cortexRegroup(new BundlesFusion((int) n));

    vector<rc_ptr<SelectBundlesFromNames> > selectNimBundles(n);
    vector<rc_ptr<SelectBundlesFromLength> > selectNBundlesFromLength(n);
    // regroup "not in mesh" bundles by gyrus name
    rc_ptr<BundlesFusion> nimRegroup( new BundlesFusion((int) n));

    // Create optional weights reader vectors
    vector<rc_ptr<FibersWeightsReader> > fibersWeightsReader(n);

    // duplicate branches for all input files
    for (i=0; i!=n; ++i) {
      // Bundles reader creation
      string fileName = fileNameIn[i];
      bundle[i].reset( new BundleReader( fileName ) );

      // Read fiber weight if file is provided
      if ( !weightsFilename.empty() ) {
        // TODO: handle multiple files given
        fibersWeightsReader[i].reset(
          new FibersWeightsReader( weightsFilename ));
        bundle[i]->addBundleListener(*fibersWeightsReader[i]);

        // Set names from mesh/label texture
        gyriFilter[i].reset(
          new SelectFiberListenerFromMesh(mesh, tex, mode, addInt, motion,""));
        fibersWeightsReader[i]->addBundleListener(*gyriFilter[i]);
      }
      else {
        // Set names from mesh/label texture
        gyriFilter[i].reset(
          new SelectFiberListenerFromMesh(mesh, tex, mode, addInt, motion,""));
        bundle[i]->addBundleListener(*gyriFilter[i]);
      }

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

    //  Set writers if weight file provided
    if ( !weightsFilename.empty() ) {

      // cortex fibers weights writer
      rc_ptr< FibersWeightsWriter > cortexWeightsWriter;
      string cortexWeightsFilename = fileNameOut + "cortexWeigths";
      cortexWeightsWriter.reset(new FibersWeightsWriter(cortexWeightsFilename));
      cortexRegroup->addBundleListener(*cortexWeightsWriter);

      // cortex bundles writer
      rc_ptr< BundleWriter > cortexWriter;
      cortexWriter.reset(new BundleWriter);
      cortexWriter->setFileString(fileNameOut);
      cortexWeightsWriter->addBundleListener(*cortexWriter);

      // NIM fibers weights writer
      rc_ptr< FibersWeightsWriter > nimWeightsWriter;
      string nimWeightsFilename = fileNameOut_notinmesh + "nimWeigths.txt";
      nimWeightsWriter.reset(new FibersWeightsWriter(nimWeightsFilename));
      nimRegroup->addBundleListener(*nimWeightsWriter);

      // NIM bundles writer
      rc_ptr< BundleWriter > nimWriter(new BundleWriter);
      nimWriter->setFileString(fileNameOut_notinmesh);
      nimWeightsWriter->addBundleListener(*nimWriter);
    }
    // Set writers if no weight file provided
    else {
      // cortex bundles writer
      rc_ptr< BundleWriter > cortexWriter;
      cortexWriter.reset(new BundleWriter);
      cortexWriter->setFileString(fileNameOut);
      cortexRegroup->addBundleListener(*cortexWriter);

      // NIM bundles writer
      rc_ptr< BundleWriter > nimWriter(new BundleWriter);
      nimWriter->setFileString(fileNameOut_notinmesh);
      nimRegroup->addBundleListener(*nimWriter);
    }

    // run all those
    for (i = 0; i < n; ++i) {
      if (verbose) cout << "process file: " << fileNameIn[i] << endl;
      bundle[i]->read();
      if (verbose) cout << "done for " << fileNameIn[i] << endl;
    }

    if (verbose) cout << "All done.\n";

  return EXIT_SUCCESS;
  }
  catch (carto::user_interruption &) {
    // Exceptions thrown by command line parser (already handled, simply exit)
  }
  catch (exception & e) {
    cerr << e.what() << endl;
  }
}
