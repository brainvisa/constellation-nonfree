#ifndef CONSTELLATION_FUSIONMESHROIS_TO_TEXTURE_H
#define CONSTELLATION_FUSIONMESHROIS_TO_TEXTURE_H

#include <aims/mesh/surface.h>
#include <aims/mesh/texture.h>
#include <aims/resampling/motion.h>

class Graph;

namespace constel {

  TimeTexture<short> *fusionMeshROIs_to_Texture_S16(
      const AimsSurfaceTriangle &inAimsMesh, const Graph &ROI_graph,
      float distthresh, std::string vertex_property = "aims_roi");

  TimeTexture<short> *fusionMeshFoldGraph_to_Texture(
      const AimsSurfaceTriangle &inAimsMesh, const Graph &ROI_graph,
      float distthresh, std::string vertex_property = "aims_roi",
      std::string label_name = "roi_label");

  Graph *texMeshAndBundles_to_BundlesGraph(
      const AimsSurfaceTriangle &inAimsMesh,
      carto::rc_ptr<TimeTexture<short> > tex, std::string namesMode,
      std::string bundlesFile_name, Motion motion,
      const std::string &BundlesNamesfile_name = std::string(),
      float fibers_filter_proportion = 0., int texture_time_step = 0);

  Graph *texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing(
      const AimsSurfaceTriangle &inAimsMesh, const TimeTexture<short> &tex,
      std::string bundlesFile_name, Motion motion,
      std::string BundlesNamesfile_name);

} // namespace constel

#endif // ifndef CONSTELLATION_FUSIONMESHROIS_TO_TEXTURE_H
