# SPDX-License-Identifier: MIT
cmake_minimum_required(VERSION 3.14)

# Subproject of efika
project(efika-apss)

# Include build configurations (compiler standards, options, etc.)
if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND EFIKA_BUILDCONFIG_FILE)
  include(${EFIKA_BUILDCONFIG_FILE})
endif()

# Include required CMake modules
include(GenerateExportHeader)
include(FetchContent)

#-------------------------------------------------------------------------------
# COMPONENT configuration
#-------------------------------------------------------------------------------
add_library(${PROJECT_NAME}     OBJECT)
add_library(${PROJECT_NAME}-sdk INTERFACE)

generate_export_header(${PROJECT_NAME}
  EXPORT_FILE_NAME efika/apss/export.h
  BASE_NAME        apss
  PREFIX_NAME      EFIKA_)

target_sources(${PROJECT_NAME}
  PRIVATE src/allpairs.c src/allpairs_pp.c src/bruteforce.c src/idxjoin.c
          src/idxjoin_pp.c #[[src/sfr.c]])

target_include_directories(${PROJECT_NAME}
  PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/include
          ${CMAKE_CURRENT_BINARY_DIR})

target_include_directories(${PROJECT_NAME}-sdk
  INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/include)

add_library(Efika::apss    ALIAS ${PROJECT_NAME})
add_library(EfikaSDK::apss ALIAS ${PROJECT_NAME}-sdk)

#-------------------------------------------------------------------------------
# Declare efika dependencies
#-------------------------------------------------------------------------------
foreach(dep core)
  FetchContent_Declare(efika-${dep}
    GIT_REPOSITORY git@github.com:jiverson002/efika-${dep}.git)

  FetchContent_MakeAvailable(efika-${dep})

  target_link_libraries(${PROJECT_NAME} PRIVATE Efika::${dep} EfikaSDK::${dep})
endforeach()

#-------------------------------------------------------------------------------
# Unit test configuration
#-------------------------------------------------------------------------------
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  include(CTest)
  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()
endif()
