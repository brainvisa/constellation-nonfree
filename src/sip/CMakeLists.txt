set( _mainSipFile "${CMAKE_CURRENT_SOURCE_DIR}/constellation.sip" )

file( GLOB _sipSources "${CMAKE_CURRENT_SOURCE_DIR}/*.sip" )

file( MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/${CONSTELLATION-NONFREE_RELATIVE_SIP_DIRECTORY}" )

BRAINVISA_ADD_SIP_PYTHON_MODULE( constelsip
  constel
  "${_mainSipFile}"
  SIP_SOURCES ${_sipSources}
  SIP_INCLUDE "${CMAKE_BINARY_DIR}/${CONSTELLATION-NONFREE_RELATIVE_SIP_DIRECTORY}"
    "${AIMS-FREE_SIP_DIRECTORY}"
  SIP_INSTALL "${CONSTELLATION-NONFREE_RELATIVE_SIP_DIRECTORY}" )
target_link_libraries( constelsip "${PYTHON_LIBRARY}" constellation pyaims )

BRAINVISA_INSTALL( TARGETS constelsip
                   DESTINATION python/constel
                   COMPONENT ${PROJECT_NAME} )

