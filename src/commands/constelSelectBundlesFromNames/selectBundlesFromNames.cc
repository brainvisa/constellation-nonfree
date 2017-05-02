#include <aims/getopt/getopt2.h>
#include <aims/resampling/linearInterpolator.h>
#include <aims/fibers/bundles_features.h>
#include <cartobase/config/verbose.h>
#include <constellation/selectBundlesFromNames.h>

using namespace std;
using namespace carto;
using namespace aims;
using namespace constel;


int main(int argc, const char **argv) {
  try {
    string bunIn;
    string bunOut;
    vector<string> namesList;
    bool verbose  = false;
    bool as_regex = false;

    AimsApplication app(
        argc, argv, "Select bundles from the input bundles file."
        "The bundle names to be selected are specified by a names list.\n");

    app.addOption(
        bunIn, "-i",
        "bundles input");
    app.addOption(
        bunOut, "-o",
        "bundles output");
    app.addOptionSeries(
        namesList, "-names",
        "list of string containing bundle names to be selected");
    app.addOption(
        verbose, "-verbose",
        "show as much information as possible", true);
    app.addOption(
        as_regex, "-r",
        "names are regular expressions", true);

    app.initialize();

    // Reader
    BundleReader bundleReader(bunIn);

    // Select bundles from names
    rc_ptr<SelectBundlesFromNames> selectBundlesFromNames;
    selectBundlesFromNames.reset(
        new SelectBundlesFromNames(namesList, verbose, as_regex));
    bundleReader.addBundleListener(*selectBundlesFromNames);

    // Writer
    rc_ptr<BundleWriter> bundleWriter;
    bundleWriter.reset(new BundleWriter());
    bundleWriter->setFileString(bunOut);
    selectBundlesFromNames->addBundleListener(*bundleWriter);

    // Running engine
    bundleReader.read();

    return EXIT_SUCCESS;
  }
  catch (user_interruption &) {}

  catch (exception & e) {
    cerr << e.what() << endl;
  }
  return EXIT_FAILURE;
}
