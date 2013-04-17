
#include <constellation/textureAndMeshTools.h>
#include <aims/mesh/surfaceOperation.h>
#include <algorithm>

using namespace std;
using namespace aims;

namespace constel
{

  template<int D, class T> std::vector<std::set<uint> >
  surfacePolygonsIndex( const AimsTimeSurface<D,T> & surf )
  {
    /*
    output:       out_vector : std::vector<std::set<uint> > of size vertex nb: out_vector[i] = index of polygons associated to vertex i in surf.polygon()
    */
    const std::vector< AimsVector<uint,D> > & poly = surf.polygon();
    uint n = poly.size();
    std::vector<std::set<uint> > polygons_index( surf.vertex().size() );

    for ( uint i=0; i<n; ++i )
    {
      for ( uint j=0; j<D; ++j )
      {
        polygons_index[poly[i][j]].insert( i );
      }
    }
    return polygons_index;
  }//surfacePolygonsIndex


  std::set<uint> surfacePolygonsIndexByVerticesGroup(const std::vector<std::set<uint> > & polygonsByVertex_Index, std::vector<uint> & vertexIndex)
  {
    // a tester
    uint n = vertexIndex.size();
    std::set<uint> polygonsByVerticesGroup_Index;
    std::set<uint>::const_iterator set_it;
    for ( uint i=0; i<n; ++i )
    {
      uint currentVertex_index = vertexIndex[i];
      for (set_it = polygonsByVertex_Index[currentVertex_index].begin(); set_it != polygonsByVertex_Index[currentVertex_index].end(); set_it++)
      {
        polygonsByVerticesGroup_Index.insert(*set_it);
      }
    }
    return polygonsByVerticesGroup_Index;
  }//surfacePolygonsIndexBySetOfVertex


  bool hasVertexLabel(const TimeTexture<short> & labeled_tex, std::set<uint> vertexIndex_set, int label)
  {
//     std::cout << "Computing hasVertexLabel method..."<< std::endl;
    const Texture<short> & labeled_tex0 = labeled_tex.begin()->second;
    std::set<uint>::const_iterator set_it, set_it_begin = vertexIndex_set.begin(), set_it_end = vertexIndex_set.end();
//     std::cout << "labeled_tex0[*set_it]:" << std::flush;
    for(set_it = set_it_begin; set_it!= set_it_end; set_it++)
    {
//       std::cout << labeled_tex0[*set_it] <<", " << std::flush;
      if (labeled_tex0[*set_it]==label)
      {
        return true;
      }
    }
//     std::cout << "Done." << std::endl;
    return false;
  }//hasVertexLabel


  bool connectedCommponent_isInside(const AimsSurfaceTriangle & aimsMesh,const TimeTexture<short> & region_tex, std::vector< std::size_t > connectedCommponent_vertexIndex)
  {
    std::vector<std::set<uint> > aimsMeshNeighbours = SurfaceManip::surfaceNeighbours( aimsMesh );
    std::size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
//     std::cout << "connectedCommponent_vertexIndex size:" << connectedCommponent_size << "." << std::endl;
    std::size_t vertex_index;
//     std::cout << "vertex_index: "<< std::flush; 
    for(std::size_t i=0;i<connectedCommponent_size;i++)
    {
      vertex_index = connectedCommponent_vertexIndex[i];
//       std::cout << vertex_index << std::flush; 
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index],0))
      {
//         std::cout << " hasVertexLabel 0, " << std::flush;
        return false;
      }
//       std::cout << ", " << std::flush;
    }
//     std::cout << "Done." << std::endl;
    return true;
    

  }//connectedCommponent_isInside


  std::size_t firstBoundaryVertex(const AimsSurfaceTriangle & aimsMesh,const TimeTexture<short> & region_tex, std::vector< std::size_t > connectedCommponent_vertexIndex)
  {
    std::vector<std::set<uint> > aimsMeshNeighbours = SurfaceManip::surfaceNeighbours( aimsMesh );
    std::size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
//     std::cout << "connectedCommponent_vertexIndex size:" << connectedCommponent_size << "." << std::endl;
    std::size_t vertex_index;
//     std::cout << "vertex_index: "<< std::flush; 
    for(std::size_t i=0;i<connectedCommponent_size;i++)
    {
      vertex_index = connectedCommponent_vertexIndex[i];
//       std::cout << vertex_index << std::flush; 
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index],0))
      {
//         std::cout << " has a neighbour with VertexLabel 0 " << std::flush;
        return vertex_index;
      }
//       std::cout << ", " << std::flush;
    }
//     std::cout << "Done." << std::endl;
    return -1;
  }//firstBoundaryVertex


  template <typename T>
      std::size_t boundaryVertexBiggestValue(const AimsSurfaceTriangle & aimsMesh,const TimeTexture<short> & region_tex, std::vector< std::size_t > connectedCommponent_vertexIndex, const TimeTexture<T> & values_tex)
  {
    const Texture<T> & values_tex0 = values_tex.begin()->second;
    std::vector<std::set<uint> > aimsMeshNeighbours = SurfaceManip::surfaceNeighbours( aimsMesh );
    std::size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
//     std::cout << "connectedCommponent_vertexIndex size:" << connectedCommponent_size << "." << std::endl;
    std::size_t vertex_index;
//     std::cout << "vertex_index: "<< std::flush;
    T max_value = 0;
    T vertex_value;
    std::size_t chosen_vertex = -1;
    for(std::size_t i=0;i<connectedCommponent_size;i++)
    {
      vertex_index = connectedCommponent_vertexIndex[i];
//       std::cout << vertex_index << std::flush;
      vertex_value = values_tex0[vertex_index];
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index],0) and vertex_value > max_value)
      {
//         std::cout << " has a neighbour with VertexLabel 0 " << std::flush;
        chosen_vertex = vertex_index;
        max_value = vertex_value;
      }
//       std::cout << ", " << std::flush;
    }
//     std::cout << "Done." << std::endl;
    return chosen_vertex;
  }//boundaryVertexBiggestValue


  bool computeIntersectionSegmentPolygon( Point3df segmentBeginPoint, Point3df segmentEndPoint, std::vector< Point3df > polygon, Point3df ** intersection_point )
  {
    //polygon.size >= 3
    if (polygon.size() < 3)
      throw runtime_error("error on the polygon dimension: polygon.size() > 3");
    //from Ray-polygon intersection, Didier Badouel, Graphics Gem I

    //Step 1: Intersecting the embedding plane:
    //Computing the normal N of the plane:
    Point3df norm = Point3df(0.0);
    norm = crossed( (polygon[1]-polygon[0]), (polygon[2]- polygon[0]) );
    int i0, i1, i2;
    std::vector<float > normal(3);
    normal[0] = norm[0];
    normal[1] = norm[1];
    normal[2] = norm[2];
    float norm_comp_max = *(max_element(normal.begin(), normal.end()));
    //
    if (norm[0] == norm_comp_max)
    {
      i0 = 0;
      i1 = 1;
      i2 = 2;
    }
    else if (norm[1] == norm_comp_max)
    {
      i0 = 1;
      i1 = 0;
      i2 = 2;
    }
    else if (norm[2] == norm_comp_max)
    {
      i0 = 2;
      i1 = 0;
      i2 = 1;
    }
    //Parametric representation of the segment (by a ray):
    //r(t) = O + D.t, such as r(0) = segmentBeginPoint and r(1) = segmentEndPoint
    Point3df O = segmentBeginPoint;
    Point3df D = segmentEndPoint - segmentBeginPoint;
    float t_max = 1.0;
    //Compute d = - Vo.N, Vo = polygon[0]
    float d = - (polygon[0].dot(norm));
    //Test if polygon and ray are parallel (norm.D = 0):
    float N_D = norm.dot(D);
    if (N_D == 0)
    {
      return false;
    }
    //t corresponding to the intersection point:
    float t = -(d + norm.dot(O))/N_D;
    if (t < 0 || t > t_max)
    {
      return false;
    }
    
    //step 2: Intersecting the polygon or not ?
    Point3df P = O + t*D;//intersection point
    float u0, u1, u2, v0, v1, v2, beta, alpha;
    u0 = P[i1] - polygon[0][i1];
    v0 = P[i2] - polygon[0][i2];
    bool inter = false;
    uint i = 2;
    do
    {
      //The polygon is viewed as (n-2) triangles.
      u1 = polygon[i-1][i1] - polygon[0][i1];
      v1 = polygon[i-1][i2] - polygon[0][i2];
      u2 = polygon[i][i1] - polygon[0][i1];
      v2 = polygon[i][i2] - polygon[0][i2];
      
      if (u1 == 0)
      {
        beta = u0/u2;
        if ((beta>=0.)&&(beta<=1.))
        {
          alpha = (v0-beta*v2)/v1;
          inter = ((alpha >=0.)&&(alpha+beta)<= 1.);
        }
      }
      else
      {
        beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
        if ((beta >= 0.)&&(beta <= 1.))
        {
          alpha = (u0 - beta*u2)/u1;
          inter = ((alpha >= 0)&&((alpha+beta) <= 1.));
        }
      }
    }  while ((!inter)&&(++i < polygon.size()));
    if (intersection_point)
    {
      if (inter)
      {
        * intersection_point = new Point3df(P);
      }
      else
      {
        * intersection_point = 0;
      }
    }
    return inter;

  }//computeIntersectionSegmentPolygon


  template std::vector<std::set<uint> >
      surfacePolygonsIndex( const AimsTimeSurface<2, Void> & surf );
  template std::vector<std::set<uint> >
      surfacePolygonsIndex( const AimsSurfaceTriangle & surf );
  template std::vector<std::set<uint> >
      surfacePolygonsIndex( const AimsSurfaceFacet & surf );

  template
      std::size_t boundaryVertexBiggestValue(const AimsSurfaceTriangle & aimsMesh,const TimeTexture<short> & region_tex, std::vector< std::size_t > connectedCommponent_vertexIndex, const TimeTexture<float> & values_tex);

  template
      std::size_t boundaryVertexBiggestValue(const AimsSurfaceTriangle & aimsMesh,const TimeTexture<short> & region_tex, std::vector< std::size_t > connectedCommponent_vertexIndex, const TimeTexture<short> & values_tex);


}//namespace constel

