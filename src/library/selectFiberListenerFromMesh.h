#ifndef CONSTELLATION_SELECTFIBERLISTENERFROMMESH_H
#define CONSTELLATION_SELECTFIBERLISTENERFROMMESH_H


#include <constellation/bundleSet.h>
#include <constellation/tildefs.h>
#include <aims/mesh/texture.h>

//------------------//
// SelectFiberListenerFromMesh //
//------------------//

namespace constel
{

/** Fibers selection by regions on a mesh (plus a label texture)
*/
class SelectFiberListenerFromMesh
  : public comist::BundleProducer, public comist::BundleListener
{
public:
  /** namesMode sould be in "NameFront_NameEnd", "Name1_Name2",
      "Name1_Name2orNotInMesh", "NameFront", "NameEnd"
  */
  SelectFiberListenerFromMesh( constel::Mesh mesh, TimeTexture<short> tex,
                               std::string namesMode, int addInt,
                               Motion motion ,
                               const std::string &bundlesNamesFileName );
  virtual ~SelectFiberListenerFromMesh();
  void setStream( std::ostream & );

protected:


  virtual void bundleStarted( const comist::BundleProducer &,
                              const comist::BundleInfo & );
  virtual void bundleTerminated( const comist::BundleProducer &,
                                 const comist::BundleInfo & );
  virtual void fiberStarted( const comist::BundleProducer &,
                             const comist::BundleInfo &,
                             const comist::FiberInfo & );
  virtual void fiberTerminated( const comist::BundleProducer &,
                                const comist::BundleInfo &,
                                const comist::FiberInfo & );
  virtual void newFiberPoint( const comist::BundleProducer &,
                              const comist::BundleInfo &,
                               const comist::FiberInfo &,
                              const comist::FiberPoint & );
  virtual void noMoreBundle( const comist::BundleProducer & );

private:

  Fiber _fiber;
  Fibers _fibers;
  BundlesSet _bundlesSet;
  std::vector<std::string> _bundles_name;
  std::string _file_name;
  std::ofstream _file_internal;
  std::ostream *_file;
  std::string _new_bundle_name;

  constel::Mesh _mesh;
  TimeTexture<short> _tex;
  std::string _namesMode;
  int _addInt;
  Motion _motion;
};

} // namespace constel

#endif

