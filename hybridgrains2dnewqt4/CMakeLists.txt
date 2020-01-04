include( CMakeSourceFiles.txt )

add_executable( hybridgrains2dnew_qt4 ${Headers} ${Sources} ${HWQ_MOC_CPP} )

set_target_properties( hybridgrains2dnew_qt4 PROPERTIES AUTOMOC TRUE )

find_package( OpenGL REQUIRED )
target_include_directories( hybridgrains2dnew_qt4 SYSTEM PRIVATE ${OPENGL_INCLUDE_DIR} )

find_package( Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED )
include( ${QT_USE_FILE} )

target_link_libraries( hybridgrains2dnew_qt4 ${QT_LIBRARIES} ${OPENGL_LIBRARIES} rigidbody2dutils hybridgrains2dnewutils hybridgrains2dnew )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets/hybridgrains2dnew ${CMAKE_CURRENT_BINARY_DIR}/assets )
