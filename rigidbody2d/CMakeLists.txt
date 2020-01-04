include( CMakeSourceFiles.txt )

add_library( rigidbody2d ${Headers} ${Sources} )

target_include_directories( rigidbody2d SYSTEM PRIVATE )

target_link_libraries( rigidbody2d scisim )
