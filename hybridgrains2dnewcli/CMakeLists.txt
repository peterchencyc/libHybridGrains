include( CMakeSourceFiles.txt )

add_executable( hybridgrains2dnew_cli ${Headers} ${Sources} )

target_link_libraries( hybridgrains2dnew_cli rigidbody2dutils hybridgrains2dnewutils hybridgrains2dnew )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets/hybridgrains2dnew ${CMAKE_CURRENT_BINARY_DIR}/assets )


