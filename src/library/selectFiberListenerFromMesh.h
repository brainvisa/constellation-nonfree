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

/** Fibers selection by regions on a mesh (plus a label texture).
    This filter just outputs regions names in an external file (or stream).
    Regions are named after their intersections with a WM mesh labelled 
    with a label texture.
*/
class SelectFiberListenerFromMesh
  : public comist::BundleProducer, public comist::BundleListener
{
public:
  /** namesMode sould be in "NameFront_NameEnd", "Name1_Name2",
      "Name1_Name2orNotInMesh", "NameFront", "NameEnd"
  */
  SelectFiberListenerFromMesh( carto::rc_ptr<constel::Mesh> mesh,
                               carto::rc_ptr<TimeTexture<short> > tex,
                               const std::string & namesMode, int addInt,
                               const Motion & motion,
                               const std::string & bundlesNamesFileName,
                               int texture_time_step = 0 );
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

  std::string fiberName( const Point3df & p1, const Point3df & p2 );

private:
  struct Private;
  Private *d;

  std::string _file_name;
  std::ofstream _file_internal;
  std::ostream *_file;

  carto::rc_ptr<constel::Mesh> _mesh;
  carto::rc_ptr<TimeTexture<short> > _tex;
  std::string _namesMode;
  int _addInt;
  Motion _motion;

  bool _fiberstarted;
  Point3df _p1;
  Point3df _p2;
  std::string _current_name;
};

} // namespace constel

#endif

