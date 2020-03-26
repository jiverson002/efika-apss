cmake_minimum_required( VERSION 3.10 )

# Subproject of efika
project( efika-impl )

# Always use '-fPIC'/'-fPIE' option.
set( CMAKE_POSITION_INDEPENDENT_CODE ON )

# Require C99 compliance from compiler
set( CMAKE_C_STANDARD 99 )
set( CMAKE_C_STANDARD_REQUIRED ON )
set( CMAKE_C_EXTENSIONS OFF )

#-------------------------------------------------------------------------------
# COMPONENT configuration
#-------------------------------------------------------------------------------
add_library( ${PROJECT_NAME} )
add_library( ${PROJECT_NAME}-friend INTERFACE )

target_sources( ${PROJECT_NAME}
  PRIVATE src/sfr.c )

target_include_directories( ${PROJECT_NAME}
  PUBLIC  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
          $<INSTALL_INTERFACE:include>
  PRIVATE src/include )

target_include_directories( ${PROJECT_NAME}-friend
  INTERFACE src/include )

add_library( Efika::impl        ALIAS ${PROJECT_NAME} )
add_library( Efika::impl-friend ALIAS ${PROJECT_NAME}-friend )

#-------------------------------------------------------------------------------
# Declare efika dependencies
#-------------------------------------------------------------------------------
include( FetchContent )

foreach( dep core )
  FetchContent_Declare( efika-${dep}
    GIT_REPOSITORY git@github.com:jiverson002/efika-${dep}.git )
  FetchContent_MakeAvailable( efika-${dep} )

  target_link_libraries( ${PROJECT_NAME} PUBLIC  Efika::${dep} )
  target_link_libraries( ${PROJECT_NAME} PRIVATE Efika::${dep}-friend )
endforeach()

#-------------------------------------------------------------------------------
# Unit test configuration
#-------------------------------------------------------------------------------
if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME )
  include( CTest )
endif()
if( CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING )
  add_subdirectory( test )
endif()
