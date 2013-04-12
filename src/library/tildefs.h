#ifndef CONSTELLATION_TILDEFS_H
#define CONSTELLATION_TILDEFS_H

#include <cathier/aims_wrap.h>
#include <cathier/kdtree.h>

namespace constellation
{

  typedef til::Mesh_N MyMesh;
  typedef til::KDTree<std::size_t, til::MeshTraits<MyMesh>::VertexCollection> MyKDTree;

}

#endif

