
#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <connectomist/fibertracking/bundleRegroup.h>
#include <constellation/selectFiberListenerFromMesh.h>
#include <constellation/selectBundlesFromNames.h>

using namespace aims;
using namespace carto;
using namespace std;
using namespace comist;
using namespace constel;



int main( int argc, const char* argv[] )
{
  try
  {
    string fileNameIn;
    string fileNameOut;
    string fileNameOut_notinmesh;
    string fileBundlesNameOut="";
    string mode="Name1_Name2orNotInMesh";
    Mesh mesh;
    TimeTexture<short> tex;
    Reader<AimsSurfaceTriangle> r;
    int addInt =0;
    Reader< TimeTexture<short> > texR;
    string mname;
    std::vector< string > namesList;
    bool as_regex = false;
    bool verbose = false;
    AffineTransformation3d motion;
    string gyrus;

    AimsApplication app( argc, argv, "Selection of tracking fibers according to a white matter mesh with label texture." );
    app.addOption( fileNameIn, "-i", "Bundle File input");
    app.addOption( fileNameOut, "-o", "Bundle File output");
    app.addOption( fileNameOut_notinmesh, "-n", "Bundle File output for \"distant fibers\"");
    app.addOption( fileBundlesNameOut, "--namesfile",
                   "BundleNames File output");
    app.addOption( r, "--mesh", "input mesh" );
    app.addOption( texR, "--tex", "labels input texture" );
    app.addOption( gyrus, "-g", "gyrus name for distant fibers filtering" );
    app.addOption( mname, "--trs", "transformation from t2 to anat", true );
    app.addOption( addInt, "--intadd", "int to add to Bundle Names", true );
    app.addOption( mode, "--mode", "mode of new bundles names : = Name1_Name2 or NameFront, or NameEnd, or NameFront_NameEnd, or Name1_Name2orNotInMesh (default)",true);
    app.addOptionSeries( namesList, "--names", "list of string containing bundle names to be selected" );
    app.addOption( as_regex, "-r", "names are regular expressions", true );
    app.addOption( verbose, "--verbose",
                   "show as much information as possible", true );
    app.initialize();
    AimsSurfaceTriangle s;
    r.read( s );
    til::Mesh1 mesh0;
    til::convert(mesh0, s);
    mesh = addNeighborsToMesh(mesh0);

    if( !mname.empty() )
    {
      Reader<AffineTransformation3d> mreader(mname);
      mreader.read(motion);
    }

    if (verbose)
      std::cout << "reading texture..." << flush;
      texR.read( tex );
    if (verbose)
    {
      std::cout << "done" << std::endl;
      std::cout << "# vertices : " << getVertices(mesh).size() << std::endl;
      std::cout << "# faces : " << getFaceIndices(mesh).size() << std::endl;
      std::cout << "Texture dim : " << tex[0].nItem() << endl;
    }
    char buff[100];
    int a = int(tex[0].item(0));
    std::sprintf(buff, "%02i", a );
    if (verbose)
    {
      std::cout << "texture test " << tex[0].item(0) << std::endl;
      std::cout << "buff test " << buff << std::endl;

      std::cout << "fileBundlesNameOut:" << fileBundlesNameOut << std::endl;
    }
    //     First Bundles reader creation
    BundleReader bundle(fileNameIn);
    SelectFiberListenerFromMesh
      fiberBundle(mesh,tex,mode,addInt,motion, "" /*fileBundlesNameOut*/);
    stringstream names_stream;
    fiberBundle.setStream( names_stream );
    if (verbose) cout << "read bundle" << fileNameIn << endl;
    bundle.addBundleListener(fiberBundle);

    //     Second Bundles reader creation
    if ( verbose ) cout << "creating second bundle reader: " << fileNameIn << "..." << endl;

//     BundleReader bundleReader2(fileNameIn);
//     if ( verbose ) cout << "  done" << endl << endl;

    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset( new BundleRegroup( "" /*fileBundlesNameOut*/) );
    bundleRegroup->setStream( names_stream );
    fiberBundle.addBundleListener(*bundleRegroup);

    rc_ptr<SelectBundlesFromNames> selectBundlesFromNames;
    selectBundlesFromNames.reset( new SelectBundlesFromNames(namesList,
      verbose, as_regex ));
    bundleRegroup->addBundleListener(*selectBundlesFromNames);

    rc_ptr< BundleWriter > bundleWriter;
    bundleWriter.reset( new BundleWriter() );
    bundleWriter->setFileString(fileNameOut);
    selectBundlesFromNames->addBundleListener( *bundleWriter );

    // -- 2nd branch

    vector<string> notinmesh_names;
    notinmesh_names.push_back( gyrus + "_notInMesh" );
    rc_ptr<SelectBundlesFromNames> selectBundlesFromNames_notinmesh;
    selectBundlesFromNames_notinmesh.reset( new SelectBundlesFromNames( notinmesh_names,
      verbose, true ));
    bundleRegroup->addBundleListener(*selectBundlesFromNames_notinmesh);

    rc_ptr< BundleWriter > bundleWriter_notinmesh;
    bundleWriter_notinmesh.reset( new BundleWriter() );
    bundleWriter_notinmesh->setFileString(fileNameOut_notinmesh);
    selectBundlesFromNames_notinmesh->addBundleListener( *bundleWriter_notinmesh );

    //
    
    if (verbose) cout << "read() and create BundlesNames file" << endl;
    bundle.read();
    if (verbose) cout << "done" << endl;

/*    if ( verbose ) cout << "reading bundle and BundlesNames file to regroup fibers in: " << fileNameOut << "..." << endl;
    bundleReader2.read();*/
//     if (verbose) cout << "done" << endl;




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
