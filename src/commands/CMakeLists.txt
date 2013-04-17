#------------------------------------------------------------------------------
# Build targets for commands
#------------------------------------------------------------------------------

set( _commands
  "AimsSparseMatrixSmoothing"
  "constelConnectivityMatrix"
  "constelFibersLengths"
  "constelSelectBundlesFromNames"
  "constelSelectFibersFromMesh"
)

foreach( _command ${_commands} )
  file( GLOB _pro "${_command}/*.pro" )
  BRAINVISA_GET_FILE_LIST_FROM_PRO( "${_pro}"
                              TARGET _target
                              SOURCES _proSources )
  set( _sources )
  foreach( _f ${_proSources} )
    set( _sources ${_sources} ${_command}/${_f} )
  endforeach()
  add_executable( ${_target} ${_sources} )
  target_link_libraries( ${_target} constellation )
  BRAINVISA_INSTALL( TARGETS ${_target}
                     DESTINATION bin
                     COMPONENT ${PROJECT_NAME} )
endforeach()

