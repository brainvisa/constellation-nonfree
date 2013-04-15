#ifndef CONSTELLATION_CONNECTIVITIES_H
#define CONSTELLATION_CONNECTIVITIES_H

#include <cathier/SparseVector.h>
#include <vector>
#include <set>

namespace constel
{

  typedef til::SparseVector<double> Connectivity;
  typedef std::vector< Connectivity > Connectivities;
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  typedef std::vector< std::pair<std::size_t, double> > QuickMap;
  typedef std::vector<QuickMap> Connection;
  typedef std::vector<Connection> BundleConnections;
  typedef boost::shared_ptr<BundleConnections> BundleConnections_shared_ptr;
  typedef std::vector<uint> ConnectionsLength;
  typedef std::vector<float> ConnectionsFloatLength;
  typedef boost::shared_ptr<ConnectionsLength> ConnectionsLength_shared_ptr;
  typedef std::vector<std::set<uint> > PolygonsByVertexIndex;
  typedef std::vector<PolygonsByVertexIndex> PolygonsByVertexIndex_collection;

}

#endif

