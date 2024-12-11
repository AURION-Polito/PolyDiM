cmake_minimum_required(VERSION 3.12)
project(PolyDiM_Project LANGUAGES CXX C VERSION 1.0)

set(ENABLE_POLYDIM_UNIT_TESTS OFF CACHE BOOL "Enable Unit Tests")
set(ENABLE_POLYDIM_EXAMPLES OFF CACHE BOOL "Enable Examples")

# Include dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim)

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM)

# Include unit tests
if(${ENABLE_POLYDIM_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/test)
endif()

# Include addim examples
if(${ENABLE_POLYDIM_EXAMPLES})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/examples)
endif()
