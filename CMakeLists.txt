cmake_minimum_required(VERSION 3.10)

#-------------------------------------------------------------------------------
# BOOTSTRAP w/ latest bootstrap/cmake scripts
#-------------------------------------------------------------------------------
include(FetchContent)
FetchContent_Declare(cmake
  GIT_REPOSITORY git@github.com:jiverson002/bootstrap-cmake.git
)
FetchContent_MakeAvailable(cmake)

#-------------------------------------------------------------------------------
# Import dependencies
#-------------------------------------------------------------------------------
set(core_GIT_REPOSITORY git@github.com:jiverson002/efika-core.git)
set(data_GIT_REPOSITORY git@github.com:jiverson002/efika-data.git)

import(core data)

#-------------------------------------------------------------------------------
# SUBPROJECT configuration
#-------------------------------------------------------------------------------
project(${BS_PROJECT_NAME}-impl LANGUAGES C CXX)

cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0079 NEW)

#-------------------------------------------------------------------------------
# LIBRARY configuration
#-------------------------------------------------------------------------------
add_library(${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_link_libraries(${BS_PROJECT_NAME} INTERFACE ${PROJECT_NAME})

add_subdirectory(src)

#-------------------------------------------------------------------------------
# Require C++17 support from compiler for tests
#-------------------------------------------------------------------------------
cxx_17()

#-------------------------------------------------------------------------------
# ...
#-------------------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(BSTest)
endif()

#-------------------------------------------------------------------------------
# Unit test configuration
#-------------------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
  add_subdirectory(test)
endif()
