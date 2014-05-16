set( _headers
  "bundleLoader.h"
  "bundleSet.h"
  "bundleTools.h"
  "connectivities.h"
  "connConversion.h"
  "connMatrix.h"
  "connMatrixTools.h"
  "lengthBundleListener.h"
  "mergeMeshAndROIs.h"
  "selectBundlesFromLength.h"
  "selectBundlesFromNames.h"
  "selectFiberListenerFromMesh.h"
  "sparseMatrixSmoothing.h"
  "textureAndMeshTools.h"
  "tildefs.h"
)

set( _sources
  "bundleTools.cc"
  "bundleTools_meshConnectionBundleListener.cc"
  "bundleTools_meshIntersectionNoSmoothingFasterBundleListener.cc"
  "bundleTools_meshIntersectionNoSmoothingBundleListener.cc"
  "bundleTools_fiberNameAccordingToMeshTextureIntersectionBundleListener.cc"
  "bundleTools_fiberNameAccordingToMeshIntersectionBundleListener.cc"
  "bundleTools_meshIntersectionBundleListener.cc"
  "connConversion.cc"
  "connMatrix.cc"
  "connMatrixTools.cc"
  "connMatrixTools2.cc"
  "lengthBundleListener.cc"
  "mergeMeshAndROIs.cc"
  "selectBundlesFromLength.cc"
  "selectBundlesFromNames.cc"
  "selectFiberListenerFromMesh.cc"
  "sparseMatrixSmoothing.cc"
  "textureAndMeshTools.cc"
)

set(CONSTELLATION_HEADERS ${_headers} PARENT_SCOPE)
BRAINVISA_COPY_AND_INSTALL_HEADERS( _headers "constellation" headersTarget )
add_library( constellation SHARED ${_sources} )
add_dependencies( constellation ${headersTarget} )
set_property( TARGET constellation PROPERTY VERSION ${${PROJECT_NAME}_VERSION} )
set( _libraries ${OLD_CONNECTOMIST-PRIVATE_LIBRARIES} )
target_link_libraries( constellation ${_libraries} )

BRAINVISA_INSTALL( TARGETS constellation
              DESTINATION lib
              COMPONENT ${PROJECT_NAME} )
