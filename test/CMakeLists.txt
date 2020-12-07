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

#-------------------------------------------------------------------------------
# GoogleTest configuration
#-------------------------------------------------------------------------------
option(BENCHMARK_ENABLE_TESTING OFF)
option(BUILD_GMOCK              OFF)

FetchContent_Declare(googletest
  GIT_REPOSITORY https://github.com/google/googletest.git)

FetchContent_MakeAvailable(googletest)

set_target_properties(gtest_main
  PROPERTIES EXCLUDE_FROM_ALL 1 EXCLUDE_FROM_DEFAULT_BUILD 1)

gtest_discover_tests(${PROJECT_NAME}-tests
  PROPERTIES TIMEOUT 3600)

#-------------------------------------------------------------------------------
# Declare efika dependencies
#-------------------------------------------------------------------------------
foreach(dep core data io)
  FetchContent_Declare(efika-${dep}
    GIT_REPOSITORY https://github.com/jiverson002/efika-${dep}.git)

  FetchContent_MakeAvailable(efika-${dep})

  target_link_libraries(${PROJECT_NAME}-tests PRIVATE Efika::${dep})
endforeach()
