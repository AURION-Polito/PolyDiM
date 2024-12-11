list(APPEND elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)

list(APPEND elliptic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)

list(APPEND elliptic_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(elliptic_PCC_2D_sources ${elliptic_PCC_2D_sources} PARENT_SCOPE)
set(elliptic_PCC_2D_headers ${elliptic_PCC_2D_headers} PARENT_SCOPE)
set(elliptic_PCC_2D_includes ${elliptic_PCC_2D_includes} PARENT_SCOPE)
