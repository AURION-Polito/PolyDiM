list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_LocalSpace.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_Inertia_LocalSpace.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_LocalSpace_Data.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_ReferenceElement.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_Ortho_LocalSpace.hpp)

list(APPEND polydim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_LocalSpace.cpp)
list(APPEND polydim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_Ortho_LocalSpace.cpp)
list(APPEND polydim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VEM_PCC_2D_Inertia_LocalSpace.cpp)

list(APPEND polydim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_src_sources ${polydim_src_sources} PARENT_SCOPE)
set(polydim_src_headers ${polydim_src_headers} PARENT_SCOPE)
set(polydim_src_includes ${polydim_src_includes} PARENT_SCOPE)
