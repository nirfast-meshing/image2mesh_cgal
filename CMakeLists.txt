# Created by the script cgal_create_cmake_script
# This is the CMake script for compiling a CGAL application.


project( image2mesh_cgal )

cmake_minimum_required(VERSION 2.6.2)
if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" VERSION_GREATER 2.6)
  if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION}" VERSION_GREATER 2.8.3)
    cmake_policy(VERSION 2.8.4)
  else()
    cmake_policy(VERSION 2.6)
  endif()
endif()

find_package(CGAL QUIET COMPONENTS Core ImageIO)

if ( CGAL_FOUND )

  include( ${CGAL_USE_FILE} )

  include( CGAL_CreateSingleSourceCGALProgram )
  #add_definitions(-O2)
  set(CMAKE_BUILD_TYPE Release)
  create_single_source_cgal_program( "image2mesh_cgal.cpp" )

  add_library(static_boost_system STATIC IMPORTED)
  set_property(TARGET static_boost_system PROPERTY IMPORTED_LOCATION /usr/local/lib/libboost_system-mt.a)
add_library(static_boost_thread STATIC IMPORTED)
  set_property(TARGET static_boost_thread PROPERTY IMPORTED_LOCATION /usr/local/lib/libboost_thread-mt.a)
add_library(static_gmpxx STATIC IMPORTED)
  set_property(TARGET static_gmpxx PROPERTY IMPORTED_LOCATION /usr/local/lib/libgmpxx.a)
add_library(static_mpfr STATIC IMPORTED)
  set_property(TARGET static_mpfr PROPERTY IMPORTED_LOCATION /usr/local/lib/libmpfr.a)
add_library(static_gmp STATIC IMPORTED)
  set_property(TARGET static_gmp PROPERTY IMPORTED_LOCATION /usr/local/lib/libgmp.a)
  target_link_libraries(image2mesh_cgal static_boost_system static_boost_thread static_gmpxx static_gmp static_mpfr )
  add_definitions(-Wl,-no-pie)
else()

    message(STATUS "This program requires the CGAL library, and will not be compiled.")

endif()

