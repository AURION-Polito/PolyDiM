cmake_minimum_required(VERSION 3.12)
project(PolyDiM_Project LANGUAGES CXX C VERSION 1.0)

set(ENABLE_GEDIM_UNIT_TESTS OFF CACHE BOOL "Enable Gedim Unit Tests")
set(ENABLE_POLYDIM_UNIT_TESTS OFF CACHE BOOL "Enable Polydim Unit Tests")
set(ENABLE_POLYDIM_EXAMPLES OFF CACHE BOOL "Enable Polydim Examples")

# Include dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim/GeDiM)

if(${ENABLE_GEDIM_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim/GeDiM/test)
endif()

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM)

# Include unit tests
if(${ENABLE_POLYDIM_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/test)
endif()

# Include examples
if(${ENABLE_POLYDIM_EXAMPLES})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/examples)
endif()

include(${CMAKE_CURRENT_SOURCE_DIR}/gedim/cmake/cppcheck.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/gedim/cmake/clang-format.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/gedim/cmake/prepend_license.cmake)
