cmake_minimum_required(VERSION 3.13)

#-------------------------------------------------------------------------------
# BOOTSTRAP
#-------------------------------------------------------------------------------
include(${CMAKE_CURRENT_LIST_DIR}/cmake/Bootstrap.cmake)

#-------------------------------------------------------------------------------
# SUBPROJECT configuration
#-------------------------------------------------------------------------------
project(${BS_PROJECT_NAME}-impl LANGUAGES C CXX)

cmake_policy(SET CMP0076 NEW)
cmake_policy(SET CMP0079 NEW)

#-------------------------------------------------------------------------------
# Require C99/C++17 support from compiler
#-------------------------------------------------------------------------------
c_99()
cxx_17()

#-------------------------------------------------------------------------------
# SUBLIBRARY configuration
#-------------------------------------------------------------------------------
add_library(${PROJECT_NAME})

target_link_libraries(${BS_PROJECT_NAME} INTERFACE ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

target_include_directories(${PROJECT_NAME} PRIVATE
  $<TARGET_PROPERTY:${BS_PROJECT_NAME}-core,INCLUDE_DIRECTORIES>
  $<TARGET_PROPERTY:${BS_PROJECT_NAME}-io,INCLUDE_DIRECTORIES>
)

add_subdirectory(src)

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
