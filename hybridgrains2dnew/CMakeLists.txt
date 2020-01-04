include( CMakeSourceFiles.txt )

add_library( hybridgrains2dnew ${Headers} ${Sources} )

target_link_libraries( hybridgrains2dnew rigidbody2d mpmgrains2d mpmgrains2dutils )
if( USE_OPENMP )
  find_package( OpenMP )
  if( NOT OPENMP_FOUND )
    message( FATAL_ERROR "Error, failed to locate OpenMP." )
  endif()
  target_compile_options(  hybridgrains2dnew PRIVATE ${OpenMP_CXX_FLAGS} )
  # target_link_libraries( scisim INTERFACE ${OpenMP_CXX_FLAGS} )
endif()
