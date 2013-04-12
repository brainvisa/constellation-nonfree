set( _headers "tildefs.h"
  "sparseMatrix.h"
  "sparseMatrixItemR.h"
  "sparseMatrixItemW.h"
)

set( _sources "sparseMatrix.cc"
  "sparseMatrixItemR.cc"
  "sparseMatrixItemW.cc"
)

set(LIB_HEADERS ${_headers} PARENT_SCOPE)
BRAINVISA_COPY_AND_INSTALL_HEADERS( _headers "constellation" headersTarget )
add_library( constellation SHARED ${_sources} )
add_dependencies( constellation ${headersTarget} )
set_property( TARGET constellation PROPERTY VERSION ${${PROJECT_NAME}_VERSION} )
set( _libraries ${OLD_CONNECTOMIST-PRIVATE_LIBRARIES} )
target_link_libraries( constellation ${_libraries} )

BRAINVISA_INSTALL( TARGETS constellation
              DESTINATION lib
              COMPONENT ${PROJECT_NAME} )
