#include <constellation/selectFiberListenerFromMesh.h>

using namespace std;
using namespace carto;
using namespace aims;
using namespace constel;

namespace constel {

  struct SelectFiberListenerFromMesh::Private
  {
    Private(const Mesh &mesh, int texture_time_step)
    : kdt(),
      tex_time_step(texture_time_step)
      {
        KDTreeVertices m = kdt_vertices( mesh );
        kdt = KDTree( m.begin(), m.end() );
      }

    KDTree kdt;
    int tex_time_step;
  };


  //-------------------------------------------------------------------------
  SelectFiberListenerFromMesh::SelectFiberListenerFromMesh(
      rc_ptr<Mesh> mesh, rc_ptr<TimeTexture<short> > tex,
      const string &namesMode, int addInt, const Motion &motion,
      const string &bundlesNamesFileName, int texture_time_step)
      : d(new Private(*mesh, texture_time_step)),
        _file_name(bundlesNamesFileName), _file(0),
        _mesh(mesh), _tex(tex),
       _namesMode(namesMode), _addInt(addInt), _motion(motion)
  {
    if (!_file_name.empty()) {
      _file_internal.open(_file_name.c_str());
      _file = &_file_internal;
      if (_file_internal.is_open())
          cout << "BundlesNames file opened: " << _file_name << endl;
    }
  }


  //-------------------------------------------------------------------------
  SelectFiberListenerFromMesh::~SelectFiberListenerFromMesh() {
    if (_file_internal.is_open()) _file_internal.close();
    delete d;
  }


  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::bundleStarted(
      const BundleProducer &, const BundleInfo & bundleInfo) {
    _current_name.clear();
    startBundle(bundleInfo);
  }

  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::bundleTerminated(
      const BundleProducer &, const BundleInfo & bundleInfo) {
    if (!_current_name.empty()) {
      BundleInfo binfo(_current_name);
      terminateBundle(binfo);
    } else {
      terminateBundle(bundleInfo);
    }
  }


  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::fiberStarted(
      const BundleProducer &, const BundleInfo & bundleInfo,
      const FiberInfo & fiberInfo) {
    _fiberstarted = false;
    startFiber(bundleInfo, fiberInfo);
  }


  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::fiberTerminated(
      const BundleProducer &, const BundleInfo & bundleInfo,
      const FiberInfo & fiberInfo) {
    string name = fiberName(_p1, _p2);
    if (_file)
      *_file << name << endl;
    if (name == "255")
      cout << "SelectFiberListenerFromMesh: fiber 255 found\n";

    // trick BundleInfo to set fiber name in
    BundleInfo bundleInfo2(bundleInfo.id(), name);
    terminateFiber(bundleInfo2, fiberInfo);
    _current_name = name;
  }


  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::newFiberPoint(
      const BundleProducer &, const BundleInfo & bundleInfo,
      const FiberInfo & fiberInfo, const FiberPoint &point) {
    if (_fiberstarted) {
      _p2 = point;
    } else {
      _fiberstarted = true;
      _p1 = point;
    }
    addFiberPoint(bundleInfo, fiberInfo, point);
  }


  //-------------------------------------------------------------------------
  string SelectFiberListenerFromMesh::fiberName(
      const Point3df &p1T2, const Point3df &p2T2) {
    Point3df p1, p2;
    p1 = _motion.transform(p1T2[0], p1T2[1], p1T2[2]);
    til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
    std::size_t A = d->kdt.find_nearest( make_pair( 0U, p1na ) ).first->first;
    std::size_t B;
    p2 = _motion.transform(p2T2[0], p2T2[1], p2T2[2]);
    til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
    B = d->kdt.find_nearest( make_pair( 0U, p2na ) ).first->first;
    stringstream fibername;
    if (_namesMode=="NameFront_NameEnd") {
      if (til::dist2(p1na, getVertices(*_mesh)[A], til::prec<float>())
          <= 25.0 &&  // fiber point is close to the mesh
          til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
          <= 25.0) {
        fibername << int((*_tex)[d->tex_time_step].item(A)) + _addInt << "_"
          << int((*_tex)[d->tex_time_step].item(B)) + _addInt;
      } else {
        fibername << "trash";
      }
    }
    if (_namesMode=="Name1_Name2") {
      float meshClosestPointMaxDistance2 = 25.0;
      if (til::dist2(p1na, getVertices(*_mesh)[A], til::prec<float>())
          <= meshClosestPointMaxDistance2 &&  // fiber point is close to the mesh, before meshClosestPointMaxDistance2 = 25.0
          til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
          <= meshClosestPointMaxDistance2) {
        int nameA = int((*_tex)[d->tex_time_step].item(A));
        int nameB = int((*_tex)[d->tex_time_step].item(B));
        if (nameA < nameB) {
          fibername << nameA + _addInt << "_" << nameB + _addInt;
        } else {
          fibername << nameB + _addInt << "_" << nameA + _addInt;
        }
      } else {
        fibername << "trash";
      }
    } else if (_namesMode=="Name1_Name2orNotInMesh") {
      //allow that an extremity is very far from cortex
      float meshClosestPointMaxDistance2 = 25.0;
      if (til::dist2(p1na, getVertices(*_mesh)[A], til::prec<float>())
          <= meshClosestPointMaxDistance2) {
        int nameA = int((*_tex)[d->tex_time_step].item(A));
        if (til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
            <= meshClosestPointMaxDistance2) {
          int nameB = int((*_tex)[d->tex_time_step].item(B));
          if (nameA < nameB) {
            fibername << nameA + _addInt << "_" << nameB + _addInt;
          } else {
            fibername << nameB + _addInt << "_" << nameA + _addInt;
          }
        } else {
          fibername << nameA + _addInt << "_notInMesh";
        }
      } else {
        if (til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
            <= meshClosestPointMaxDistance2) {
          int nameB = int((*_tex)[d->tex_time_step].item(B));
          fibername << nameB + _addInt << "_notInMesh";
        } else {
          fibername << "trash";
        }
      }
    } else if(_namesMode=="NameFront") {
      if (til::dist2(p1na, getVertices(*_mesh)[A], til::prec<float>())
          <= 25.0 &&  // fiber point is close to the mesh
          til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
          <= 25.0) {
        fibername << int((*_tex)[d->tex_time_step].item(A)) + _addInt;
      } else {
        fibername << "trash";
      }
    } else if(_namesMode=="NameEnd") {
      if (til::dist2(p1na, getVertices(*_mesh)[A], til::prec<float>())
          <= 25.0 &&  // fiber point is close to the mesh
          til::dist2(p2na, getVertices(*_mesh)[B], til::prec<float>())
          <= 25.0) {
        fibername << int((*_tex)[d->tex_time_step].item(B)) + _addInt;
      } else {
        fibername << "trash";
      }
    }
    return fibername.str();
  }

  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::noMoreBundle(const BundleProducer &) {
    BundleProducer::noMoreBundle();
  }

  //-------------------------------------------------------------------------
  void SelectFiberListenerFromMesh::setStream(ostream &stm) {
    _file = &stm;
  }

} //namespace constel

