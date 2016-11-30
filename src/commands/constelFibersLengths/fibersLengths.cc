#include <aims/getopt/getopt2.h>
#include <constellation/lengthBundleListener.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace constel;

int main(int argc, const char* argv[]) {
  try {
    string bundleFilename;
    string outputFibersLengthFilename;
    bool verbose = false;

    AimsApplication app(
        argc, argv, "Output Fibers Lengths in a text file");

    app.addOption(
        bundleFilename, "-i",
        "input bundles" );
    app.addOption(
        outputFibersLengthFilename, "-o",
        "output fibers_length text filename");
    app.addOption(
        verbose, "-verbose",
        "show as much information as possible", true);

    app.initialize();
    
    BundleReader bundleReader(bundleFilename);
    LengthBundleListener lengthBundleListener(outputFibersLengthFilename);
    bundleReader.addBundleListener(lengthBundleListener);
    bundleReader.read();

    return EXIT_SUCCESS;
  }
  catch (carto::user_interruption &) {
    // Exceptions thrown by command line parser (already handled, simply exit)
  }
  catch (exception & e) {
    cerr << e.what() << endl;
  }
  return EXIT_FAILURE;
}

