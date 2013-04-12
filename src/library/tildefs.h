#ifndef CONSTELLATION_TILDEFS_H
#define CONSTELLATION_TILDEFS_H

#include <cathier/aims_wrap.h>
#include <cathier/kdtree.h>

namespace constel
{

  typedef til::Mesh_N Mesh;
  typedef til::KDTree<std::size_t, til::MeshTraits<Mesh>::VertexCollection> KDTree;

}

#endif

