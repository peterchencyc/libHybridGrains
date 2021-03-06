include( CMakeSourceFiles.txt )

# Track the Git commit within the source
include( GetGitRevisionDescription )
get_git_head_revision( GIT_REFSPEC GIT_SHA1 )
list( APPEND Sources "${CMAKE_CURRENT_BINARY_DIR}/CompileDefinitions.cpp" )

# Create a file with build information
configure_file( "${CMAKE_CURRENT_SOURCE_DIR}/CompileDefinitions.cpp.in" "${CMAKE_CURRENT_BINARY_DIR}/CompileDefinitions.cpp" @ONLY )

if( CMAKE_Fortran_COMPILER )
  set_source_files_properties( CompileDefinitions.h CompileDefinitions.cpp PROPERTIES COMPILE_DEFINITIONS FORTRAN_FOUND )
endif()

add_library( scisim ${Headers} ${Sources} ${FORTRAN_SOURCES} )

# Add a reliable indicator of the current compiler
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  target_compile_definitions( scisim PUBLIC CMAKE_DETECTED_GNU_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  target_compile_definitions( scisim PUBLIC CMAKE_DETECTED_INTEL_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  target_compile_definitions( scisim PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  target_compile_definitions( scisim PUBLIC CMAKE_DETECTED_CLANG_COMPILER )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()

# Eigen3 is header only and required for callers of the scisim library
find_package( Eigen3 3.2.7 REQUIRED )
target_include_directories( scisim SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR} )
# Disable Eigen's built in multithreading as we will use OpenMP around Eigen
target_compile_definitions( scisim PUBLIC EIGEN_DONT_PARALLELIZE )
if( CMAKE_BUILD_TYPE STREQUAL "Debug" )
  # EIGEN_DONT_ALIGN also enabled EIGEN_DONT_VECTORIZE
  target_compile_definitions( scisim PUBLIC EIGEN_DONT_ALIGN EIGEN_INITIALIZE_MATRICES_BY_NAN )
endif()

# MKL is an external library and is used internal to scisim, but required for linking
if( USE_MKL )
  find_package( MKL REQUIRED )
  target_include_directories( scisim SYSTEM PUBLIC ${MKL_INCLUDE_DIRS} )
  target_compile_definitions( scisim PUBLIC MKL_FOUND )
  target_link_libraries( scisim INTERFACE ${MKL_LIBRARIES} )
  # MKL requires pthreads on Linux with GCC
  if( NOT APPLE )
    if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
      find_package( Threads )
      target_link_libraries( scisim INTERFACE ${CMAKE_THREAD_LIBS_INIT} )
    endif()
  endif()
endif()

# HDF5 is an external library and is used internal to scisim
if( USE_HDF5 )
  find_package( HDF5 REQUIRED )
  if( NOT ${HDF5_FOUND} )
    message( FATAL_ERROR "HDF5 not found." )
  endif()
  target_include_directories( scisim SYSTEM PUBLIC ${HDF5_INCLUDE_DIRS} )
  target_compile_definitions( scisim PUBLIC USE_HDF5 )
  target_link_libraries( scisim INTERFACE ${HDF5_LIBRARIES} )
endif()

# Embedded Python and the Interpreter are external libraries and required for callers of scisim
if( USE_PYTHON )
  find_package( PythonLibs 2 REQUIRED )
  target_include_directories( scisim SYSTEM PUBLIC ${PYTHON_INCLUDE_DIRS} )
  find_package( NumPy REQUIRED )
  target_include_directories( scisim SYSTEM PUBLIC ${NUMPY_INCLUDE_DIRS} )
  target_compile_definitions( scisim PUBLIC USE_PYTHON )
  target_link_libraries( scisim INTERFACE ${PYTHON_LIBRARIES} )
endif()

# OpenMP is only used in the core scisim library but required when linking to scisim
if( USE_OPENMP )
  find_package( OpenMP )
  if( NOT OPENMP_FOUND )
    message( FATAL_ERROR "Error, failed to locate OpenMP." )
  endif()
  target_compile_options( scisim PRIVATE ${OpenMP_CXX_FLAGS} )
  target_link_libraries( scisim INTERFACE ${OpenMP_CXX_FLAGS} )
endif()
