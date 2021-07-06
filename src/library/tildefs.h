#ifndef CONSTELLATION_TILDEFS_H
#define CONSTELLATION_TILDEFS_H

#include <cathier/aims_wrap.h>
#include <kdtree++/kdtree.hpp>

namespace constel {
  typedef til::Mesh_N Mesh;

  template <typename PointType>
  struct Bracket_accessor_PointIndex
  {
    typedef float result_type;
    typedef std::pair<uint, PointType> _Val;

    result_type
    operator()(_Val const& V, size_t const N) const
    {
      return V.second[N];
    }
  };

  typedef std::vector< std::pair< uint, til::MeshTraits<Mesh>::Vertex > >
    KDTreeVertices;

  /// get vertices vector with index for each vertex
  inline
  KDTreeVertices kdt_vertices( const Mesh & mesh )
  {
    std::vector<std::pair< uint, til::MeshTraits<Mesh>::Vertex > > vert;
    vert.reserve( getVertices(mesh).size() );
    til::MeshTraits<Mesh>::VertexCollection::const_iterator
      iv, ev = getVertices(mesh).end();
    uint index = 0;
    for( iv=getVertices(mesh).begin(); iv!=ev; ++iv, ++index )
      vert.push_back( std::make_pair( index, *iv ) );
    return vert;
  }

  typedef KDTree::KDTree<3, std::pair<uint,
    til::MeshTraits<Mesh>::Vertex>,
    Bracket_accessor_PointIndex<til::MeshTraits<Mesh>::Vertex> >
    KDTree;
}
#endif

