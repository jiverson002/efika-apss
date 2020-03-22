# Build tests executable
add_executable(${PROJECT_NAME}-tests
  test_interface.cc
)

# ...
target_link_libraries(${PROJECT_NAME}-tests
  ${BS_PROJECT_NAME}
)

# Discover tests.
project_add_tests(${PROJECT_NAME}-tests
  PROPERTIES TIMEOUT 20
)
