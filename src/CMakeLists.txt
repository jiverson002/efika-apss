#-------------------------------------------------------------------------------
# LIBRARY configuration
#-------------------------------------------------------------------------------
target_sources(${PROJECT_NAME} PRIVATE
  sfr.c
)

target_include_directories(${PROJECT_NAME} PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/include
)
