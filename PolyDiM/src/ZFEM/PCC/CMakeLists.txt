add_subdirectory(2D)

list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/ZFEM_PCC_Utilities.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/ZFEM_PCC_PerformanceAnalysis.hpp)


list(APPEND polydim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_src_sources ${polydim_src_sources} PARENT_SCOPE)
set(polydim_src_headers ${polydim_src_headers} PARENT_SCOPE)
set(polydim_src_includes ${polydim_src_includes} PARENT_SCOPE)
