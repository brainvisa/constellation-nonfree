#include <aims/getopt/getopt2.h>
#include <aims/resampling/linearInterpolator.h>
#include <connectomist/features/bundles_features.h>
#include <cartobase/config/verbose.h>

#include <constellation/selectBundlesFromNames.h>


using namespace std;
using namespace carto;
using namespace aims;
using namespace comist;
using namespace constel;


int main( int argc, const char **argv )
{
  try
  {
    AimsApplication app( argc, argv, "Select bundles from the input bundles file."
        "The bundle names to be selected are specified by a names list.\n" );
    string bunIn;
    string bunOut;
    std::vector< string > namesList;
    bool verbose = false;
    bool as_regex = false;

    app.addOption( bunIn, "-i", "bundles input" );
    app.addOption( bunOut, "-o", "bundles output" );
    app.addOptionSeries( namesList, "-names", "list of string containing bundle names to be selected" );
    app.addOption( verbose, "-verbose", "show as much information as possible", true );
    app.addOption( as_regex, "-r", "names are regular expressions", true );
    app.initialize();
    if ( verbose ) cout<<endl;

    if ( verbose ) cout << "Creating bundle reader: " << bunIn << "..." << flush;
    BundleReader bundleReader(bunIn);
    if ( verbose ) cout << "  done" << endl;

    rc_ptr<SelectBundlesFromNames> selectBundlesFromNames;
    selectBundlesFromNames.reset( new SelectBundlesFromNames(namesList,
      verbose, as_regex ));
    bundleReader.addBundleListener(*selectBundlesFromNames);

    rc_ptr< BundleWriter > bundleWriter;
    bundleWriter.reset( new BundleWriter() );
    bundleWriter->setFileString(bunOut);
    selectBundlesFromNames->addBundleListener( *bundleWriter );

    // Running engine
    if ( verbose ) cout << "Reading bundle:..." << flush;
    bundleReader.read();
    if ( verbose ) cout << "  done" << endl;

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
