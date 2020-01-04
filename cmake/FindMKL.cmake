# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC         :   use static linking
#   MKL_MULTI_THREADED :   use multi-threading
#   MKL_SDL            :   Single Dynamic Library interface
#
# This module defines the following variables:
#
#   MKL_FOUND            : true if MKL is found
#   MKL_INCLUDE_DIRS     : where to find mkl.h, etc.
#   MKL_LIBRARIES        : the library to link against.

# TODO: Automate the generation of this file based on
#   https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

include( FindPackageHandleStandardArgs )

unset( MKL_FOUND CACHE )
unset( MKL_INCLUDE_DIRS CACHE )
unset( MKL_LIBRARIES CACHE )

set( INTEL_ROOT "/opt/intel" )
set( MKL_ROOT ${INTEL_ROOT}/mkl )

# Build the path to the libraries on Windows, OS X, and Unix
if( WIN32 )
  message( ERROR "WIN32 not yet tested for MKL, please email the author." )
elseif( APPLE )
  set( MKL_LIBRARY_PATH ${MKL_ROOT}/lib )
else() # UNIX
  set( MKL_LIBRARY_PATH ${MKL_ROOT}/lib/intel64 )
endif()

# Find the include directory
find_path( MKL_INCLUDE_DIRS mkl.h PATHS ${MKL_ROOT}/include )

# Find include directory
#  There is no include folder under Linux
if( WIN32 )
  find_path( INTEL_INCLUDE_DIR omp.h PATHS ${INTEL_ROOT}/include )
  set( MKL_INCLUDE_DIRS ${MKL_INCLUDE_DIRS} ${INTEL_INCLUDE_DIR} )
endif()

# Find libraries

# Set the suffix for static or dynamic
set( MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES} )
if( WIN32 )
  if( MKL_STATIC )
    set( CMAKE_FIND_LIBRARY_SUFFIXES .lib )
  else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES _dll.lib )
  endif()
elseif( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
  if( MKL_STATIC )
    set( CMAKE_FIND_LIBRARY_SUFFIXES .a )
  else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES .dylib )
  endif()
else() # Unix
  if( MKL_STATIC )
    set( CMAKE_FIND_LIBRARY_SUFFIXES .a )
  else()
    set( CMAKE_FIND_LIBRARY_SUFFIXES .so )
  endif()
endif()


# MKL is composed by four layers: Interface, Threading, Computational and RTL

if( MKL_SDL )
  find_library( MKL_LIBRARIES mkl_rt PATHS ${MKL_LIBRARY_PATH} )
  set( MKL_MINIMAL_LIBRARY ${MKL_LIBRARIES} )
else()
  # Set the interface layer name
  if( WIN32 )
    set( MKL_INTERFACE_LIBNAME mkl_intel_c )
  else()
    set( MKL_INTERFACE_LIBNAME mkl_intel_lp64 )
  endif()
  # Find the interface layer
  find_library( MKL_INTERFACE_LIBRARY ${MKL_INTERFACE_LIBNAME} PATHS ${MKL_LIBRARY_PATH} )

  # Set the threading layer name
  if( MKL_MULTI_THREADED )
    set( MKL_THREADING_LIBNAME mkl_intel_thread )
  else()
    set( MKL_THREADING_LIBNAME mkl_sequential )
  endif()
  # Find the threading layer
  find_library( MKL_THREADING_LIBRARY ${MKL_THREADING_LIBNAME} PATHS ${MKL_LIBRARY_PATH} )

  # Find the computational layer
  find_library( MKL_CORE_LIBRARY mkl_core PATHS ${MKL_LIBRARY_PATH} )
  #find_library( MKL_FFT_LIBRARY mkl_cdft_core PATHS ${MKL_LIBRARY_PATH} )
  #find_library( MKL_SCALAPACK_LIBRARY mkl_scalapack_core PATHS ${MKL_LIBRARY_PATH} )

  if( MKL_MULTI_THREADED )
    # Set the RTL layer name
    if( WIN32 )
      set( MKL_RTL_LIBNAME libiomp5md )
      set( MKL_RTL_PATH ${INTEL_ROOT} )
    elseif( APPLE )
      set( MKL_RTL_LIBNAME iomp5 )
      set( MKL_RTL_PATH ${INTEL_ROOT}/lib )
    else() # UNIX
      set( MKL_RTL_LIBNAME iomp5 )
      set( MKL_RTL_PATH ${INTEL_ROOT}/lib/intel64 )
    endif()

    # Find the RTL layer
    find_library( MKL_RTL_LIBRARY ${MKL_RTL_LIBNAME} PATHS ${MKL_RTL_PATH} )
  endif()

  set( MKL_LIBRARIES ${MKL_INTERFACE_LIBRARY} ${MKL_THREADING_LIBRARY} ${MKL_CORE_LIBRARY} ${MKL_RTL_LIBRARY} CACHE FILEPATH "Path to a library." )
endif()

set( CMAKE_FIND_LIBRARY_SUFFIXES ${MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES} )

find_package_handle_standard_args( MKL DEFAULT_MSG MKL_INCLUDE_DIRS MKL_LIBRARIES )

if( UNIX AND NOT APPLE )
  if( MKL_STATIC )
    SET( MKL_LIBRARIES "-Wl,--start-group" ${MKL_LIBRARIES} )
    SET( MKL_LIBRARIES ${MKL_LIBRARIES} "-Wl,--end-group" )
  endif()
endif()

unset( MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES )
unset( MKL_CORE_LIBRARY CACHE )
unset( MKL_INTERFACE_LIBRARY CACHE )
unset( MKL_THREADING_LIBRARY CACHE )
unset( MKL_RTL_LIBRARY CACHE )
mark_as_advanced( MKL_INCLUDE_DIRS )
mark_as_advanced( MKL_LIBRARIES )
