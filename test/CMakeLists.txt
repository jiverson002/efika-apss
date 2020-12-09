# Include required CMake modules
include(GoogleTest)

#-------------------------------------------------------------------------------
# EXECUTABLE configuration
#-------------------------------------------------------------------------------
add_executable(${PROJECT_NAME}-tests)

target_sources(${PROJECT_NAME}-tests
  PRIVATE test_interface.cc)

target_link_libraries(${PROJECT_NAME}-tests
  PRIVATE gtest ${PROJECT_NAME})

if(BUILD_SHARED_LIBS)
  # XXX: On Windows platform with Visual Studio, the executable and Google Test
  # dll are put in different directories and at runtime, the loader cannot find
  # the dll. One solution is to force CMake to build dll's and exe's in
  # the same directory so at runtime the dll's can be found.
  #
  # https://stackoverflow.com/a/47260387
  set_target_properties(${PROJECT_NAME}-tests
    PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

  # This should be defined whenever Google Test is built as a shared library.
  #
  # https://github.com/google/googletest/blob/master/googletest/README.md
  target_compile_definitions(${PROJECT_NAME}-tests
    PRIVATE GTEST_LINKED_AS_SHARED_LIBRARY=1)
endif()

#-------------------------------------------------------------------------------
# GoogleTest configuration
#-------------------------------------------------------------------------------
option(BUILD_GMOCK OFF)

# Prevent GoogleTest from overriding our compiler/linker options when building
# with Visual Studio
#
# https://github.com/google/googletest/blob/master/googletest/README.md
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git)

FetchContent_MakeAvailable(googletest)

set_target_properties(gtest_main
  PROPERTIES EXCLUDE_FROM_ALL 1)

gtest_discover_tests(${PROJECT_NAME}-tests
  PROPERTIES TIMEOUT 3600)

#-------------------------------------------------------------------------------
# Declare efika dependencies
#-------------------------------------------------------------------------------
foreach(dep core data io)
  FetchContent_Declare(efika-${dep}
    GIT_REPOSITORY https://github.com/jiverson002/efika-${dep}.git)

  FetchContent_MakeAvailable(efika-${dep})

  target_link_libraries(${PROJECT_NAME}-tests
    PRIVATE Efika::${dep} $<$<NOT:$<STREQUAL:${dep},data>>:EfikaSDK::${dep}>)
endforeach()
