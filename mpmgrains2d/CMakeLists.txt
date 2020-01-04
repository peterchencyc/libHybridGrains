include( CMakeSourceFiles.txt )

# Track the Git commit within the source
include( GetGitRevisionDescription )
get_git_head_revision( GIT_REFSPEC GIT_SHA1 )

# Track the build time
string( TIMESTAMP BUILD_DATE_TIME "%Y-%m-%d %H:%M:%S" )

# Create a file with build information
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h.in" "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.h" @ONLY )

# if( CMAKE_Fortran_COMPILER )
#   set_source_files_properties( CompileDefinitions.h CompileDefinitions.cpp PROPERTIES COMPILE_DEFINITIONS FORTRAN_FOUND )
# endif()

# if( USE_HDF5 )
#   set( Headers "${Headers};HDF5File.h" )
#   set( Sources "${Sources};HDF5File.cpp" )
# endif()

add_library( mpmgrains2d ${Headers} ${Sources} )

# Add a reliable indicator of the current compiler
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  target_compile_definitions( mpmgrains2d PUBLIC CMAKE_DETECTED_GNU_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  target_compile_definitions( mpmgrains2d PUBLIC CMAKE_DETECTED_INTEL_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  target_compile_definitions( mpmgrains2d PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  target_compile_definitions( mpmgrains2d PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()

# Eigen3 is header only and required for callers
find_package( Eigen3 3.2.9 REQUIRED )
target_include_directories( mpmgrains2d SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )
# Disable Eigen's built in multithreading as we will use OpenMP around Eigen
target_compile_definitions( mpmgrains2d PUBLIC EIGEN_DONT_PARALLELIZE )
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  # EIGEN_DONT_ALIGN also enables EIGEN_DONT_VECTORIZE
  target_compile_definitions( mpmgrains2d PUBLIC EIGEN_DONT_ALIGN EIGEN_INITIALIZE_MATRICES_BY_NAN )
endif()

if( USE_HDF5 )
  find_package( HDF5 REQUIRED C )
  if( NOT ${HDF5_FOUND} )
    message( FATAL_ERROR "HDF5 not found." )
  endif()
  target_include_directories( mpmgrains2d SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS} )
  target_compile_definitions( mpmgrains2d PUBLIC USE_HDF5 )
  target_link_libraries( mpmgrains2d INTERFACE ${HDF5_C_LIBRARIES} )
endif()

target_link_libraries( mpmgrains2d LINK_PRIVATE scisim )

# Embedded Python and the Interpreter are external libraries and required for callers of scisim
# if( USE_PYTHON )
#   find_package( PythonLibs 2 REQUIRED )
#   target_include_directories( scisim SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS} )
#   find_package( NumPy REQUIRED )
#   target_include_directories( scisim SYSTEM PUBLIC ${NUMPY_INCLUDE_DIRS} )
#   target_compile_definitions( scisim PUBLIC USE_PYTHON )
#   target_link_libraries( scisim INTERFACE ${PYTHON_LIBRARIES} )
# endif()

# OpenMP is only used in the core scisim library but required when linking to scisim
if( USE_OPENMP )
  find_package( OpenMP )
  if( NOT OPENMP_FOUND )
    message( FATAL_ERROR "Error, failed to locate OpenMP." )
  endif()
  target_compile_options( mpmgrains2d PRIVATE ${OpenMP_CXX_FLAGS} )
  target_link_libraries(mpmgrains2d PUBLIC "${OpenMP_CXX_FLAGS}")
  # target_link_libraries( scisim INTERFACE ${OpenMP_CXX_FLAGS} )
endif()
