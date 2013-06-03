#include <aims/getopt/getopt2.h>
#include <aims/io/reader.h>
#include <connectomist/fibertracking/bundleRegroup.h>
#include <constellation/selectFiberListenerFromMesh.h>

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
    string fileBundlesNameOut="";
    string mode="Name1_Name2";
    rc_ptr<Mesh> mesh;
    rc_ptr<TimeTexture<short> > tex;
    Reader<AimsSurfaceTriangle> r;
    int addInt =0;
    Reader< TimeTexture<short> > texR;
    string mname;
    bool verbose = false;
    AffineTransformation3d motion;
    AimsApplication app( argc, argv, "Selection of tracking fibers according to a white matter mesh with label texture." );
    app.addOption( fileNameIn, "-i", "Bundle File input");
    app.addOption( fileNameOut, "-o", "Bundle File output");
    app.addOption( fileBundlesNameOut, "-namesfile",
                   "BundleNames File output");
    app.addOption( r, "-mesh", "input mesh" );
    app.addOption( texR, "-tex", "labels input texture" );
    app.addOption( mname, "-trs", "transformation from t2 to anat", true );
    app.addOption( addInt, "-intadd", "int to add to Bundle Names", true );
    app.addOption( mode, "-mode", "mode of new bundles names : = Name1_Name2 (default) or NameFront, or NameEnd, or NameFront_NameEnd, or Name1_Name2orNotInMesh",true);
    app.addOption( verbose, "-verbose",
                   "show as much information as possible", true );
    app.initialize();
    rc_ptr<AimsSurfaceTriangle> s;
    s.reset( r.read() );
    til::Mesh1 mesh0;
    til::convert(mesh0, *s);
    mesh.reset( new Mesh );
    *mesh = addNeighborsToMesh(mesh0);

    if( !mname.empty() )
    {
      Reader<AffineTransformation3d> mreader(mname);
      mreader.read(motion);
    }

    if (verbose)
      std::cout << "reading texture..." << flush;
      tex.reset( texR.read() );
    if (verbose)
    {
      std::cout << "done" << std::endl;
      std::cout << "# vertices : " << getVertices(*mesh).size() << std::endl;
      std::cout << "# faces : " << getFaceIndices(*mesh).size() << std::endl;
      std::cout << "Texture dim : " << (*tex)[0].nItem() << endl;
    }
    if (verbose)
      std::cout << "fileBundlesNameOut:" << fileBundlesNameOut << std::endl;

    //     First Bundles reader creation
    BundleReader bundle(fileNameIn);
    SelectFiberListenerFromMesh
      fiberBundle(mesh,tex,mode,addInt,motion,fileBundlesNameOut);
    if (verbose) cout << "read bundle" << fileNameIn << endl;
    bundle.addBundleListener(fiberBundle);
    if (verbose) cout << "read() and create BundlesNames file" << endl;
    bundle.read();
    if (verbose) cout << "done" << endl;

    //     Second Bundles reader creation
    if ( verbose ) cout << "creating second bundle reader: " << fileNameIn << "..." << endl;

    BundleReader bundleReader2(fileNameIn);
    if ( verbose ) cout << "  done" << endl << endl;

    rc_ptr< BundleRegroup > bundleRegroup;
    bundleRegroup.reset( new BundleRegroup(fileBundlesNameOut) );
    bundleReader2.addBundleListener(*bundleRegroup);

    rc_ptr< BundleWriter > bundleWriter;
    bundleWriter.reset( new BundleWriter() );
    bundleWriter->setFileString(fileNameOut);
    bundleRegroup->addBundleListener( *bundleWriter );

    if ( verbose ) cout << "reading bundle and BundlesNames file to regroup fibers in: " << fileNameOut << "..." << endl;
    bundleReader2.read();
    if (verbose) cout << "done" << endl;




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
