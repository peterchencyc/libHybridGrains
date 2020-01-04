include( CMakeSourceFiles.txt )

add_library( rigidbody2dutils ${Headers} ${Sources} )

find_package( RapidXML REQUIRED )
target_include_directories( rigidbody2dutils SYSTEM PRIVATE ${RAPIDXML_INCLUDE_DIR} )

if( USE_CAIRO )
  find_package(Cairo REQUIRED)
  if( NOT CAIRO_FOUND )
    message( FATAL_ERROR "Error, failed to locate the Cairo library." )
  endif()
  target_link_libraries( rigidbody2dutils ${CAIRO_LIBRARIES} )
  target_include_directories( rigidbody2dutils SYSTEM PUBLIC ${CAIRO_INCLUDE_DIRS} )
  target_compile_definitions( rigidbody2dutils PUBLIC USE_CAIRO )
endif()

target_link_libraries( rigidbody2dutils rigidbody2d )
