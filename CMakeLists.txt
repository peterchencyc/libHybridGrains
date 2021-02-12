# 3.2.1 to allow set_property( TARGET target PROPERTY CXX_STANDARD 11 ) and similar
cmake_minimum_required( VERSION 3.2.1 )

# Check that the Eigen is installed
execute_process( COMMAND ./scripts_include/eigen_found.sh RESULT_VARIABLE EIGEN_NOT_FOUND WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
if( EIGEN_NOT_FOUND )
  message( FATAL_ERROR "Failed to locate Eigen or checksum is incorrect. Please run get_dependencies.sh." )
endif()

# Check that RapidXml is installed
execute_process( COMMAND ./scripts_include/rapidxml_found.sh RESULT_VARIABLE RAPIDXML_NOT_FOUND WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} )
if( RAPIDXML_NOT_FOUND )
  message( FATAL_ERROR "Failed to locate RapidXml or checksum is incorrect. Please run get_dependencies.sh." )
endif()

project( libHybridGrains C CXX )

################################################################################
### Configure and check the build setting options

# If the user provided a CMAKE_BUILD_TYPE, create a cache variable
if( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE )
  # Otherwise default to Release
else( NOT CMAKE_BUILD_TYPE STREQUAL "" )
  set( CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." FORCE )
endif( NOT CMAKE_BUILD_TYPE STREQUAL "" )

# Enforce the use of a valid build type
set( VALID_BUILD_TYPES None Debug Release RelWithDebInfo MinSizeRel Coverage )
if( NOT ";${VALID_BUILD_TYPES};" MATCHES ";${CMAKE_BUILD_TYPE};" )
  message( FATAL_ERROR "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel Coverage." )
endif()


################################################################################
### Set flags used across all builds

# Disable compiler language extensions
set( CMAKE_CXX_EXTENSIONS OFF )
# Enable C++14
set( CMAKE_CXX_STANDARD 14 )
set( CMAKE_CXX_STANDARD_REQUIRED ON )

# Append an extra space if flags are already set
if( CMAKE_CXX_FLAGS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} " )
endif( CMAKE_CXX_FLAGS )

# Allow the user to optionally enable compiler sanitizers
set( SANITIZER "none" CACHE STRING "Enable compiler sanitizer modes, options are: none, address, undefined" )

# Enable compiler specific flags
if( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  include( CMakeListsGNUCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel" )
  include( CMakeListsIntelCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  include( CMakeListsClangCXXFlags.txt )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
  # For now, same options to Clang and AppleClang. May change in the future.
  include( CMakeListsClangCXXFlags.txt )
else()
  message( FATAL_ERROR "Unsupported C++ compiler ${CMAKE_CXX_COMPILER_ID}. Please email the maintainer." )
endif()


################################################################################
### Set paths, etc

# Create a default location for custom CMake modules
set( CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake )

# Allow projects to reference source relative to the root of the project
include_directories( ${CMAKE_CURRENT_SOURCE_DIR} )


################################################################################
### Add options to allow the user to configure the build

option( STRICT_BUILD "Enable aggressive warnings and treat warnings as errors" OFF )

option( USE_QT4 "Build Qt4 front ends" OFF )

option( USE_OPENMP "OpenMP-based parallelism" OFF )

option( USE_HDF5 "Serialization support with HDF5" OFF )

option( USE_PYTHON "Embedded Python language scripting" OFF )

option( USE_TIMERS "Enable timers" OFF )

option( USE_CAIRO "Enable suppot for rendering from the 2D command line front-ends with Cairo (https://www.cairographics.org)." OFF )

option( USE_MKL "Use Intel's MKL" OFF )
option( MKL_MULTI_THREADED "Use multithreaded MKL" OFF )
option( MKL_STATIC "Statically link against MKL" OFF )

option( USE_GPROF "Use gprof" OFF )
option( USE_VTUNE "Use vtune" OFF )

###############################################################################
### Check for installed Python modules

execute_process( COMMAND python -c "import h5py" RESULT_VARIABLE H5PY_MISSING ERROR_QUIET )
execute_process( COMMAND python -c "import numpy" RESULT_VARIABLE NUMPY_MISSING ERROR_QUIET )


###############################################################################
### Add libraries, plugins
if( USE_TIMERS )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_TIMERS" )
endif()

if( USE_OPENMP )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DOPENMP_ENABLED" )
endif()

if( USE_GPROF )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")
endif()

if( USE_VTUNE )
  set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g" )
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -g")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -g")
endif()

# Allow targets to generate tests
enable_testing()


# Core SCISim library
add_subdirectory( scisim )

# Core two dimensional rigid body library
add_subdirectory( rigidbody2d )

# Utilities shared by the command-line and Qt4 interface
add_subdirectory( rigidbody2dutils )

# Command line interface for two dimensional rigid body simulations
add_subdirectory( rigidbody2dcli )

# Qt4 interface for two dimensional rigid body simulations
if( USE_QT4 )
  add_subdirectory( rigidbody2dqt4 )
endif()

add_subdirectory( mpmgrains2d )
add_subdirectory( mpmgrains2dutils )

# Core hybrid grains simulation library
add_subdirectory( hybridgrains2dnew )
add_subdirectory( hybridgrains2dnewutils )
# Qt4 interface for the hybrid grains simulation
if( USE_QT4 )
  add_subdirectory( hybridgrains2dnewqt4 )
endif()

# Command line interface for hybrid grains simulation
add_subdirectory( hybridgrains2dnewcli )