#ifndef CONSTELLATION_CONNECTIVITIES_H
#define CONSTELLATION_CONNECTIVITIES_H

#include <vector>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>

namespace constel
{

  class Connectivity : public std::map<size_t, double>
  {
  public:
    Connectivity( size_t size ) : _size( size ) {}
    size_t size() const { return _size; }

    float get( size_t i ) const
    {
      const_iterator it = find( i );
      if( it == end() )
        return 0.;
      return it->second;
    }

    Connectivity & operator += ( const Connectivity & c )
    {
      Connectivity::const_iterator i, e = c.end();
      iterator j = find( i->first );
      if( j == end() )
        (*this)[ i->first ] = i->second;
      else
        j->second += i->second;
    }

    double dot( const Connectivity & c ) const
    {
      double sum = 0.;
      if( size() > c.size() )
      {
        Connectivity::const_iterator i, j, e = c.end(), je = end();
        for( i=c.begin(); i!=e; ++i )
        {
          j = find( i->first );
          if( j != je )
            sum += i->second * i->second;
        }
      }
      else
      {
        Connectivity::const_iterator i, j, e = end(), je = c.end();
        for( i=begin(); i!=e; ++i )
        {
          j = c.find( i->first );
          if( j != je )
            sum += i->second * i->second;
        }
      }
      return sum;
    }

  private:
    size_t _size;
  };

  typedef std::vector<Connectivity> Connectivities;
  typedef std::vector<std::vector<std::size_t> > CNeighborhoods;
  typedef std::vector<std::pair<std::size_t, double> > QuickMap;
  typedef std::map<std::size_t, float> DistMap;
  typedef std::vector<QuickMap> Connection;
  typedef std::vector<Connection> BundleConnections;
  typedef boost::shared_ptr<BundleConnections> BundleConnections_shared_ptr;
  typedef std::vector<unsigned> ConnectionsLength;
  typedef std::vector<float> ConnectionsFloatLength;
  typedef boost::shared_ptr<ConnectionsLength> ConnectionsLength_shared_ptr;
  typedef std::vector<std::set<unsigned> > PolygonsByVertexIndex;
  typedef std::vector<PolygonsByVertexIndex> PolygonsByVertexIndex_collection;

}
#endif

