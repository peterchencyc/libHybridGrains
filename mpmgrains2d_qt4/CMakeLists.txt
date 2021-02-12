include( CMakeSourceFiles.txt )

add_executable( mpmgrains2d_qt4 ${Headers} ${Sources} ${HWQ_MOC_CPP} )

set_target_properties( mpmgrains2d_qt4 PROPERTIES AUTOMOC TRUE )

find_package( OpenGL REQUIRED )
target_include_directories( mpmgrains2d_qt4 SYSTEM PRIVATE ${OPENGL_INCLUDE_DIR} )

find_package( Qt4 COMPONENTS QtCore QtGui QtOpenGL REQUIRED )
include( ${QT_USE_FILE} )

target_link_libraries( mpmgrains2d_qt4 ${QT_LIBRARIES} ${OPENGL_LIBRARIES} mpmgrains2d mpmgrains2dutils scisim )

execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_SOURCE_DIR}/assets/mpmgrains2d ${CMAKE_CURRENT_BINARY_DIR}/assets )
