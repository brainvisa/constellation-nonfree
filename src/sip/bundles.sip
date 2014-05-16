%ModuleHeaderCode
#include <constel/mergeMeshAndROIs.h>
%End

namespace constel {

    Graph *
        texMeshAndBundles_to_BundlesGraph( const AimsSurfaceTriangle & inAimsMesh, rc_ptr_TimeTexture_S16 tex, std::string namesMode, std::string bundlesFile_name, Motion motion, const std::string & BundlesNamesfile_name = std::string() ) /Factory/;

}; // namespace constel