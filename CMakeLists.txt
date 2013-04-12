cmake_minimum_required( VERSION 2.6 )
BRAINVISA_PROJECT()

find_package( brainvisa-cmake REQUIRED )

BRAINVISA_FIND_PACKAGE( old_connectomist-private REQUIRED )
BRAINVISA_FIND_PACKAGE( cortical_surface-private REQUIRED )

find_package( python )
if( PYTHON_FOUND )
  find_package( SIP )
  find_package( Numpy )
endif()

include( "${CORTICAL_SURFACE-PRIVATE_USE_FILE}" )
include( "${OLD_CONNECTOMIST-PRIVATE_USE_FILE}" )

# add_subdirectory( src/library )
# add_subdirectory( src/commands )
# if( SIP_FOUND )
#   add_subdirectory( src/sip )
# endif()

BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/python"
                                 ${PROJECT_NAME} )
BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/bin"
                          bin
                          ${PROJECT_NAME} )

