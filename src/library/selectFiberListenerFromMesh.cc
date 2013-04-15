#include <constellation/selectFiberListenerFromMesh.h>

using namespace std;
using namespace carto;
using namespace aims;
using namespace comist;
using namespace constel;

//------------------//
// SelectFiberListenerFromMesh //
//------------------//

namespace constel
{

//-----------------------------------------------------------------------------
SelectFiberListenerFromMesh::SelectFiberListenerFromMesh( Mesh mesh, TimeTexture<short> tex, string namesMode, int addInt, Motion motion, const string &bundlesNamesFileName ) : _mesh(mesh), _tex(tex), _namesMode(namesMode), _addInt(addInt), _motion(motion), _file_name(bundlesNamesFileName), _file( 0 )
{
  if( !_file_name.empty() )
  {
    _file_internal.open(_file_name.c_str() );
    _file = &_file_internal;
    if (_file_internal.is_open()) cout << "BundlesNames file opened: " << _file_name << endl;
  }
}


//-----------------------------------------------------------------------------
SelectFiberListenerFromMesh::~SelectFiberListenerFromMesh()
{
  if( _file_internal.is_open() )
    _file_internal.close();
}


//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::bundleStarted( const BundleProducer &,
                                   const BundleInfo & )
{
  _fibers.clear();
}

//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::bundleTerminated( const BundleProducer &,
                                      const BundleInfo & )
{
  _bundlesSet.push_back(_fibers);
}


//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::fiberStarted( const BundleProducer &,
                                  const BundleInfo &,
                                  const FiberInfo & )
{
  _fiber.clear();
}


//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::fiberTerminated( const BundleProducer &, 
                                     const BundleInfo &,
                                     const FiberInfo & )
{
  _fibers.push_back(_fiber);
}


//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::newFiberPoint( const BundleProducer &, 
                                   const BundleInfo &,
                                   const FiberInfo &, 
                                   const FiberPoint &point )
{
  _fiber.push_back(point);
}


//-----------------------------------------------------------------------------
void SelectFiberListenerFromMesh::noMoreBundle( const BundleProducer & )
{
//   std::cout << "Generating kdtree" << std::endl;
  KDTree kdt(getVertices(_mesh));
  makeKDTree(getVertices(_mesh), kdt);
  til::Find_closest< double, KDTree > fc(kdt);

//   output file writing
  for (BundlesSet::iterator bun = _bundlesSet.begin(); bun!=_bundlesSet.end(); ++bun)
  {
    for (Fibers::iterator fib=(*bun).begin(); fib!=(*bun).end(); ++fib)
    {
      Point3df p1T2, p2T2, p1, p2;
      p1T2 = (*fib)[0]; // fiber front
      p1 = _motion.transform(p1T2[0], p1T2[1], p1T2[2]);
      til::numeric_array<float, 3> p1na(p1[0], p1[1], p1[2]);
      std::size_t A = fc(p1na);
      std::size_t B;
      p2T2 = (*fib)[(*fib).size()-1];  // fiber end
      p2 = _motion.transform(p2T2[0], p2T2[1], p2T2[2]);
      til::numeric_array<float, 3> p2na(p2[0], p2[1], p2[2]);
      B = fc(p2na);
      if(_namesMode=="NameFront_NameEnd")
      {
        if (til::dist2(p1na, getVertices(_mesh)[A], til::prec<float>()) <= 25.0 &&  // fiber point is close to the mesh
                    til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= 25.0)
        {
          if( _file )
            *_file << int(_tex[0].item(A)) + _addInt << "_" << int(_tex[0].item(B)) + _addInt << endl;
        }
        else
        {
          if( _file )
            *_file << "trash" << endl;
        }
      }
      if(_namesMode=="Name1_Name2")
      {
        float meshClosestPointMaxDistance2 = 25.0;
//         meshClosestPointMaxDistance2 = 100.0;
        if (til::dist2(p1na, getVertices(_mesh)[A], til::prec<float>()) <= meshClosestPointMaxDistance2 &&  // fiber point is close to the mesh, before meshClosestPointMaxDistance2 = 25.0
                    til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= meshClosestPointMaxDistance2)
        {
          int nameA = int(_tex[0].item(A));
          int nameB = int(_tex[0].item(B));
          if( _file )
          {
            if (nameA < nameB)
              *_file << nameA + _addInt << "_" << nameB + _addInt << endl;
            else
            {
              *_file << nameB + _addInt << "_" << nameA + _addInt << endl;
            }
          }
        }
        else
        {
          if( _file )
            *_file << "trash" << endl;
        }
      }
      else if (_namesMode=="Name1_Name2orNotInMesh")
      {
        //allow that an extremity is very far from cortex
        float meshClosestPointMaxDistance2 = 25.0;
//         meshClosestPointMaxDistance2 = 100.0;
        if (til::dist2(p1na, getVertices(_mesh)[A], til::prec<float>()) <= meshClosestPointMaxDistance2 )
        {
          int nameA = int(_tex[0].item(A));
          if (til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= meshClosestPointMaxDistance2)
          {
            int nameB = int(_tex[0].item(B));
//             std::cout << nameB << endl;
            if( _file )
            {
              if (nameA < nameB)
              *_file << nameA + _addInt << "_" << nameB + _addInt << endl;
              else
              {
                *_file << nameB + _addInt << "_" << nameA + _addInt << endl;
              }
            }
          }
          else
          {
            if( _file )
              *_file << nameA + _addInt << "_notInMesh" << endl;
          }
        }
        else
        {
          if (til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= meshClosestPointMaxDistance2)
          {
            int nameB = int(_tex[0].item(B));
            if( _file )
              *_file << nameB + _addInt << "_notInMesh" << endl;
          }
          else
          {
            if( _file )
              *_file << "trash" << endl;
          }
        }
      }
      else if(_namesMode=="NameFront")
      {
        if (til::dist2(p1na, getVertices(_mesh)[A], til::prec<float>()) <= 25.0 &&  // fiber point is close to the mesh
                    til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= 25.0)
        {
          if( _file )
            *_file << int(_tex[0].item(A)) + _addInt << endl;
        }
        else
        {
          if( _file )
            *_file << "trash" << endl;
        }
      }
      else if(_namesMode=="NameEnd")
      {
        if (til::dist2(p1na, getVertices(_mesh)[A], til::prec<float>()) <= 25.0 &&  // fiber point is close to the mesh
                    til::dist2(p2na, getVertices(_mesh)[B], til::prec<float>()) <= 25.0)
        {
          if( _file )
            *_file << int(_tex[0].item(B)) + _addInt << endl;
        }
        else
        {
          if( _file )
            *_file << "trash" << endl;
        }
      }
    }
  }
//   cout<<endl<<"fin de l'écriture!!"<<endl<<flush;
}


void SelectFiberListenerFromMesh::setStream( ostream & stm )
{
  _file = &stm;
}

} //namespace constel
