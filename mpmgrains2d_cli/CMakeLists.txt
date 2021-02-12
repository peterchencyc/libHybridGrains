include( CMakeSourceFiles.txt )

add_executable( mpmgrains2d_cli ${Headers} ${Sources} )

target_link_libraries( mpmgrains2d_cli mpmgrains2dutils mpmgrains2d )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets/mpmgrains2d ${CMAKE_CURRENT_BINARY_DIR}/assets )

#if( USE_HDF5 )
#  add_test( mpm2d_serialization_00 assets/shell_scripts/execute_serialization_test.sh #assets/initial_translation_test.xml 0.000260 20 08 76800 )
#  add_test( mpm2d_serialization_01 assets/shell_scripts/execute_serialization_test.sh #assets/initial_rotation_test.xml 0.000260 20 08 76800 )
#endif()
