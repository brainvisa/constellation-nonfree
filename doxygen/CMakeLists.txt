IF ( DOXYGEN_FOUND )
  set(component_name "constellation")
  SET( DOXYFILE_HTML_HEADER "\"${CMAKE_CURRENT_SOURCE_DIR}/head_${component_name}.html\"" )
  set( _doxygenInput ${CONSTELLATION_HEADERS} )
  set( DOXYFILE_PREDEFINED "${CONSTELLATION-PRIVATE_DEFINITIONS}")
  set(aims-free_version "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set(aims_version "${aims-free_VERSION_MAJOR}.${aims-free_VERSION_MINOR}")
  set( DOXYFILE_TAGFILES
    "${CMAKE_BINARY_DIR}/share/doc/connectomist-${aims_version}/doxygen/connectomist.tag=../../connectomist-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/aimsalgo-${aims_version}/doxygen/aimsalgo.tag=../../aimsalgo-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/aimsdata-${aims_version}/doxygen/aimsdata.tag=../../aimsdata-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/graph-${aims_version}/doxygen/graph.tag=../../graph-${aims_version}/doxygen ${CMAKE_BINARY_DIR}/share/doc/cartobase-${aims_version}/doxygen/cartobase.tag=../../cartobase-${aims_version}/doxygen${CMAKE_BINARY_DIR}/share/doc/cartodata-${aims_version}/doxygen/cartodata.tag=../../cartodata-${aims_version}/doxygen")
  BRAINVISA_GENERATE_DOXYGEN_DOC(
    _doxygenInput  head_${component_name}.html ${component_name}.png ${component_name}.png
    INPUT_PREFIX "${CMAKE_BINARY_DIR}/include/${component_name}"
    COMPONENT "${component_name}" )
  add_dependencies( ${component_name}-doxygen connectomist-doxygen aimsalgo-doxygen aimsdata-doxygen graph-doxygen cartobase-doxygen cartodata-doxygen )
ENDIF ( DOXYGEN_FOUND )
