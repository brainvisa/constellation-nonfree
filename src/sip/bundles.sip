%ModuleHeaderCode
#include <constellation/mergeMeshAndROIs.h>
%End

namespace constel {

    Graph *
        texMeshAndBundles_to_BundlesGraph( const AimsSurfaceTriangle & inAimsMesh, rc_ptr_TimeTexture_S16 tex, std::string namesMode, std::string bundlesFile_name, Motion motion, const std::string & BundlesNamesfile_name = std::string(), float = 0. ) /Factory/;

    Graph *
        texMeshAndBundles_to_BundlesGraph_WithIntersectionComputing( const AimsSurfaceTriangle & inAimsMesh, const TimeTexture_S16 & tex,  std::string bundlesFile_name, Motion motion, std::string BundlesNamesfile_name) /Factory/;

    TimeTexture_S16 *
        fusionMeshROIs_to_Texture_S16(const AimsSurfaceTriangle & inAimsMesh, const Graph & ROI_graph, float distthresh, std::string vertex_property) /Factory/;

    TimeTexture_S16 *
        fusionMeshFoldGraph_to_Texture(const AimsSurfaceTriangle & inAimsMesh, const Graph & ROI_graph, float distthresh, std::string vertex_property, std::string label_name) /Factory/;

}; // namespace constel