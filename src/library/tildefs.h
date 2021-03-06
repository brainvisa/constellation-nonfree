#ifndef CONSTELLATION_TILDEFS_H
#define CONSTELLATION_TILDEFS_H

#include <kdtree++/kdtree.hpp>

/* This file, named "tildefs" is a remain of the former dependecy on the
   obsolete aims-til library. It was used mainly for distance maps and KDTree.
   These algorithms also exist in aims, so we have finally replaced them and do
   not depend on "til" any longer.

*/

namespace constel
{

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

  typedef std::vector< std::pair< uint, Point3df > > KDTreeVertices;

  /// get vertices vector with index for each vertex
  inline
  KDTreeVertices kdt_vertices( const AimsSurfaceTriangle & mesh )
  {
    std::vector<std::pair< uint, Point3df > > vert;
    vert.reserve( mesh.vertex().size() );
    std::vector<Point3df>::const_iterator iv, ev = mesh.vertex().end();
    uint index = 0;
    for( iv=mesh.vertex().begin(); iv!=ev; ++iv, ++index )
      vert.push_back( std::make_pair( index, Point3df( *iv ) ) );
    return vert;
  }

  typedef KDTree::KDTree<3, std::pair<uint, Point3df>,
    Bracket_accessor_PointIndex<Point3df> >
    KDTree;

  inline
  float dist2( const Point3df & p1, const Point3df & p2 )
  {
    return (p2 - p1).norm2();
  }

}
#endif

