cmake_minimum_required(VERSION 3.12)
project(PolyDiM_Project LANGUAGES CXX C VERSION 1.0)

set(ENABLE_UNIT_TESTS OFF CACHE BOOL "Enable Unit Tests")

# Include dependencies
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/gedim)

# Include library
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM)

# Include unit tests
if(${ENABLE_UNIT_TESTS})
  add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/test)
endif()
