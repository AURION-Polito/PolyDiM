list(APPEND elliptic_PCC_1D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND elliptic_PCC_1D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND elliptic_PCC_1D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND elliptic_PCC_1D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND elliptic_PCC_1D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND elliptic_PCC_1D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND elliptic_PCC_1D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(elliptic_PCC_1D_sources ${elliptic_PCC_1D_sources} PARENT_SCOPE)
set(elliptic_PCC_1D_headers ${elliptic_PCC_1D_headers} PARENT_SCOPE)
set(elliptic_PCC_1D_includes ${elliptic_PCC_1D_includes} PARENT_SCOPE)
