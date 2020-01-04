include( CMakeSourceFiles.txt )

add_library( mpmgrains2dutils ${Headers} ${Sources} )
# if( ENABLE_IWYU )
#   set_property( TARGET mpmgrains2dutils PROPERTY CXX_INCLUDE_WHAT_YOU_USE ${iwyu_path} )
# endif()

find_package( RapidXML REQUIRED )
target_include_directories( mpmgrains2dutils SYSTEM PRIVATE ${RAPIDXML_INCLUDE_DIR} )

target_link_libraries( mpmgrains2dutils mpmgrains2d )
