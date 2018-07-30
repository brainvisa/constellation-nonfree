cmake_minimum_required( VERSION 2.6 )
BRAINVISA_PROJECT()

find_package( brainvisa-cmake REQUIRED )

BRAINVISA_FIND_PACKAGE( aims-free REQUIRED )
math( EXPR result "${aims-free_VERSION_MINOR}+1" )
BRAINVISA_DEPENDENCY( RUN DEPENDS "aims-free" RUN ">= ${aims-free_VERSION};<< ${aims-free_VERSION_MAJOR}.${result}" )
BRAINVISA_DEPENDENCY( DEV DEPENDS "aims-free" DEV ">= ${aims-free_VERSION};<< ${aims-free_VERSION_MAJOR}.${result}" )

find_package( python )
if( PYTHON_FOUND )
  find_package( SIP )
  find_package( Numpy )
  if( UMPY_FOUND AND SIP_FOUND )
    BRAINVISA_DEPENDENCY( RUN DEPENDS python RUN )
    BRAINVISA_DEPENDENCY( DEV DEPENDS python DEV )
    BRAINVISA_DEPENDENCY( RUN DEPENDS python-numpy RUN )
    BRAINVISA_DEPENDENCY( RUN DEPENDS python-sip4 RUN )
    BRAINVISA_DEPENDENCY( DEV DEPENDS python-sip4 DEV )
  endif()
endif()
find_package( Doxygen )
find_package( anatomist-free )

include( "${AIMS-FREE_USE_FILE}" )

set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall" )

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
  math( EXPR result "${aims-free_VERSION_MINOR}+1" )
  BRAINVISA_DEPENDENCY( RUN RECOMMENDS "anatomist-gpl" RUN ">= ${anatomist-gpl_VERSION};<< ${anatomist-gpl_VERSION_MAJOR}.${result}" )
  BRAINVISA_COPY_PYTHON_DIRECTORY(
    "${CMAKE_CURRENT_SOURCE_DIR}/share/anatomist"
    ${PROJECT_NAME}
    "share/anatomist-${anatomist-free_VERSION_MAJOR}.${anatomist-free_VERSION_MINOR}"
  )
endif()

BRAINVISA_COPY_DIRECTORY( "${CMAKE_CURRENT_SOURCE_DIR}/bin"
                          bin
                          ${PROJECT_NAME} )

BRAINVISA_CREATE_CMAKE_CONFIG_FILES()
