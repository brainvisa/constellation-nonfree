#ifndef CONSTELLATION_TEXTUREANDMESHTOOLS_H
#define CONSTELLATION_TEXTUREANDMESHTOOLS_H

#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>

namespace constel {

  template<int D, class T> std::vector<std::set<unsigned> >
    surfacePolygonsIndex(const AimsTimeSurface<D,T> & surf);

  std::set<unsigned> surfacePolygonsIndexByVerticesGroup(
      const std::vector<std::set<unsigned> > &polygonsByVertex_Index,
      std::vector<unsigned> &vertexIndex);

  bool hasVertexLabel(
      const TimeTexture<short> &labeled_tex,
      std::set<unsigned> vertexIndex_set, int label);

  bool connectedCommponent_isInside(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      std::vector< std::size_t > connectedCommponent_vertexIndex);

  std::size_t firstBoundaryVertex(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      std::vector< std::size_t > connectedCommponent_vertexIndex);

  template <typename T> std::size_t boundaryVertexBiggestValue(
      const AimsSurfaceTriangle &aimsMesh,
      const TimeTexture<short> &region_tex,
      std::vector<std::size_t> connectedCommponent_vertexIndex,
      const TimeTexture<T> &values_tex);

  bool computeIntersectionSegmentPolygon(
      Point3df segmentBeginPoint, Point3df segmentEndPoint,
      std::vector<Point3df> polygon, Point3df **intersection_point = 0);

} // namespace constel

#endif // ifndef CONSTELLATION_TEXTUREANDMESHTOOLS_H
