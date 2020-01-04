include( CMakeSourceFiles.txt )

add_library( hybridgrains2dnewutils ${Headers} ${Sources} )

find_package( RapidXML REQUIRED )
target_include_directories( hybridgrains2dnewutils SYSTEM PRIVATE ${RAPIDXML_INCLUDE_DIR} )

target_link_libraries( hybridgrains2dnewutils hybridgrains2dnew )
