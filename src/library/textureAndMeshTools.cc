
#include <constellation/textureAndMeshTools.h>
#include <aims/mesh/surfaceOperation.h>
#include <algorithm>

using namespace std;
using namespace aims;

namespace constel {

  //------------------------
  //  surfacePolygonsIndex
  //------------------------
  template<int D, class T> vector<set<unsigned> > surfacePolygonsIndex(
      const AimsTimeSurface<D,T> &surf) {
    /*
    out_vector : vector<set<unsigned> > of size vertex nb: out_vector[i]
                 = index of polygons associated to vertex i in surf.polygon()
    */
    const vector< AimsVector<unsigned,D> > &poly = surf.polygon();
    unsigned n = poly.size();
    vector<set<unsigned> > polygons_index(surf.vertex().size());

    for (unsigned i = 0; i < n; ++i) {
      for (unsigned j = 0; j < D; ++j) {
        polygons_index[poly[i][j]].insert(i);
      }
    }
    return polygons_index;
  }

  //---------------------------------------
  //  surfacePolygonsIndexByVerticesGroup
  //---------------------------------------
  set<unsigned> surfacePolygonsIndexByVerticesGroup(
      const vector<set<unsigned> > &polygonsByVertex_Index,
      vector<unsigned> &vertexIndex) {
    // a tester
    unsigned n = vertexIndex.size();
    set<unsigned> polygonsByVerticesGroup_Index;
    set<unsigned>::const_iterator set_it;
    for (unsigned i=0; i<n; ++i) {
      unsigned currentVertex_index = vertexIndex[i];
      for (set_it = polygonsByVertex_Index[currentVertex_index].begin();
           set_it != polygonsByVertex_Index[currentVertex_index].end();
           set_it++) {
        polygonsByVerticesGroup_Index.insert(*set_it);
      }
    }
    return polygonsByVerticesGroup_Index;
  }

  //------------------
  //  hasVertexLabel
  //------------------
  bool hasVertexLabel(
      const TimeTexture<short> &labeled_tex, set<unsigned> vertexIndex_set,
      int label) {
    const Texture<short> &labeled_tex0 = labeled_tex.begin()->second;
    set<unsigned>::const_iterator set_it, set_it_begin
      = vertexIndex_set.begin(), set_it_end = vertexIndex_set.end();
    for(set_it = set_it_begin; set_it!= set_it_end; set_it++) {
      if (labeled_tex0[*set_it] == label) {
        return true;
      }
    }
    return false;
  }

  //--------------------------------
  //  connectedCommponent_isInside
  //--------------------------------
  bool connectedCommponent_isInside(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      vector<size_t> connectedCommponent_vertexIndex) {
    vector<set<unsigned> > aimsMeshNeighbours
      = SurfaceManip::surfaceNeighbours(aimsMesh);
    size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
    size_t vertex_index;
    for (size_t i=0;i<connectedCommponent_size;i++) {
      vertex_index = connectedCommponent_vertexIndex[i];
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index],0)) {
        return false;
      }
    }
    return true;
    

  }

  //-----------------------
  //  firstBoundaryVertex
  //-----------------------
  size_t firstBoundaryVertex(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      vector<size_t> connectedCommponent_vertexIndex) {
    vector<set<unsigned> > aimsMeshNeighbours
      = SurfaceManip::surfaceNeighbours(aimsMesh);
    size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
    size_t vertex_index;
    for(size_t i=0;i<connectedCommponent_size;i++) {
      vertex_index = connectedCommponent_vertexIndex[i];
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index],0)) {
        return vertex_index;
      }
    }
    return -1;
  }

  //------------------------------
  //  boundaryVertexBiggestValue
  //------------------------------
  template <typename T> size_t boundaryVertexBiggestValue(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      vector<size_t> connectedCommponent_vertexIndex,
      const TimeTexture<T> &values_tex) {
    const Texture<T> &values_tex0 = values_tex.begin()->second;
    vector<set<unsigned> > aimsMeshNeighbours
      = SurfaceManip::surfaceNeighbours(aimsMesh);
    size_t connectedCommponent_size = connectedCommponent_vertexIndex.size();
    size_t vertex_index;
    T max_value = 0;
    T vertex_value;
    size_t chosen_vertex = -1;
    for (size_t i = 0; i < connectedCommponent_size; i++) {
      vertex_index = connectedCommponent_vertexIndex[i];
      vertex_value = values_tex0[vertex_index];
      if (hasVertexLabel(region_tex, aimsMeshNeighbours[vertex_index], 0)
          and vertex_value > max_value) {
        chosen_vertex = vertex_index;
        max_value = vertex_value;
      }
    }
    return chosen_vertex;
  }

  //-------------------------------------
  //  computeIntersectionSegmentPolygon
  //-------------------------------------
  bool computeIntersectionSegmentPolygon(
      Point3df segmentBeginPoint, Point3df segmentEndPoint,
      vector<Point3df> polygon, Point3df **intersection_point) {
    //polygon.size >= 3
    if (polygon.size() < 3)
      throw runtime_error(
          "error on the polygon dimension: polygon.size() > 3");
    //from Ray-polygon intersection, Didier Badouel, Graphics Gem I

    //Step 1: Intersecting the embedding plane:
    //Computing the normal N of the plane:
    Point3df norm = Point3df(0.0);
    norm = crossed((polygon[1]-polygon[0]), (polygon[2] - polygon[0]));
    int /* i0 = 0, */i1 = 1, i2 = 2;
    vector<float > normal(3);
    normal[0] = norm[0];
    normal[1] = norm[1];
    normal[2] = norm[2];
    float norm_comp_max = *(max_element(normal.begin(), normal.end()));
    //
    if (norm[0] == norm_comp_max) {
      //i0 = 0;
      i1 = 1;
      i2 = 2;
    } else if (norm[1] == norm_comp_max) {
      //i0 = 1;
      i1 = 0;
      i2 = 2;
    } else if (norm[2] == norm_comp_max) {
      //i0 = 2;
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
    if (N_D == 0) {
      return false;
    }
    //t corresponding to the intersection point:
    float t = -(d + norm.dot(O))/N_D;
    if (t < 0 || t > t_max) {
      return false;
    }
    
    //step 2: Intersecting the polygon or not ?
    Point3df P = O + t*D; //intersection point
    float u0, u1, u2, v0, v1, v2, beta, alpha;
    u0 = P[i1] - polygon[0][i1];
    v0 = P[i2] - polygon[0][i2];
    bool inter = false;
    unsigned i = 2;
    do {
      //The polygon is viewed as (n-2) triangles.
      u1 = polygon[i-1][i1] - polygon[0][i1];
      v1 = polygon[i-1][i2] - polygon[0][i2];
      u2 = polygon[i][i1] - polygon[0][i1];
      v2 = polygon[i][i2] - polygon[0][i2];
      
      if (u1 == 0) {
        beta = u0/u2;
        if ((beta>=0.) && (beta<=1.)) {
          alpha = (v0-beta*v2)/v1;
          inter = ((alpha >=0.) && (alpha+beta)<= 1.);
        }
      } else {
        beta = (v0*u1 - u0*v1)/(v2*u1 - u2*v1);
        if ((beta >= 0.) && (beta <= 1.)) {
          alpha = (u0 - beta*u2)/u1;
          inter = ((alpha >= 0) && ((alpha+beta) <= 1.));
        }
      }
    }  while ((!inter) && (++i < polygon.size()));
    if (intersection_point) {
      if (inter) {
        *intersection_point = new Point3df(P);
      } else {
        *intersection_point = 0;
      }
    }
    return inter;
  }

  //------------------------
  //  surfacePolygonsIndex
  //------------------------
  template vector<set<unsigned> > surfacePolygonsIndex(
      const AimsTimeSurface<2, Void> &surf);

  //------------------------
  //  surfacePolygonsIndex
  //------------------------
  template vector<set<unsigned> > surfacePolygonsIndex(
      const AimsSurfaceTriangle &surf);

  //------------------------
  //  surfacePolygonsIndex
  //------------------------
  template vector<set<unsigned> > surfacePolygonsIndex(
      const AimsSurfaceFacet &surf);

  //------------------------------
  //  boundaryVertexBiggestValue
  //------------------------------
  template size_t boundaryVertexBiggestValue(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      vector<size_t> connectedCommponent_vertexIndex,
      const TimeTexture<float> &values_tex);

  //------------------------------
  //  boundaryVertexBiggestValue
  //------------------------------
  template size_t boundaryVertexBiggestValue(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      vector<size_t> connectedCommponent_vertexIndex,
      const TimeTexture<short> &values_tex);


}//namespace constel

