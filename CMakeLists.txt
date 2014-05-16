cmake_minimum_required( VERSION 2.6 )
BRAINVISA_PROJECT()

find_package( brainvisa-cmake REQUIRED )

BRAINVISA_FIND_PACKAGE( old_connectomist-private REQUIRED )

find_package( python )
if( PYTHON_FOUND )
  find_package( SIP )
  find_package( Numpy )
endif()
find_package( Doxygen )
find_package( anatomist-free )

include( "${OLD_CONNECTOMIST-PRIVATE_USE_FILE}" )

add_subdirectory( src/library )
add_subdirectory( src/commands )
if( SIP_FOUND AND PYTHON_FOUND AND NUMPY_FOUND AND SIP_FOUND )
  set( CONSTELLATION-PRIVATE_RELATIVE_SIP_DIRECTORY
    "share/constellation-${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}/sip" )
  add_subdirectory( src/sip )
endif()
add_subdirectory( doxygen )

BRAINVISA_COPY_PYTHON_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/python"
                                 ${PROJECT_NAME} )
if( anatomist-free_FOUND )
  BRAINVISA_COPY_PYTHON_DIRECTORY(
    "${CMAKE_CURRENT_SOURCE_DIR}/share/anatomist"
    ${PROJECT_NAME}
    "share/anatomist-${anatomist-free_VERSION_MAJOR}.${anatomist-free_VERSION_MINOR}"
  )
endif()

BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/bin"
                          bin
                          ${PROJECT_NAME} )
BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/R"
                          R
                          ${PROJECT_NAME} )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()
