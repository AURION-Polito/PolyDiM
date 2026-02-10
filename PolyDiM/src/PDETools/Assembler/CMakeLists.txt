list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/Assembler_Utilities.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler_PCC_2D_functions.hpp)

list(APPEND polydim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_src_sources ${polydim_src_sources} PARENT_SCOPE)
set(polydim_src_headers ${polydim_src_headers} PARENT_SCOPE)
set(polydim_src_includes ${polydim_src_includes} PARENT_SCOPE)
