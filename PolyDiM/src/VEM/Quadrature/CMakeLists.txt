list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_Quadrature_2D.hpp)
list(APPEND polydim_src_headers ${CMAKE_CURRENT_SOURCE_DIR}/VEM_Quadrature_3D.hpp)

list(APPEND polydim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VEM_Quadrature_2D.cpp)
list(APPEND polydim_src_sources ${CMAKE_CURRENT_SOURCE_DIR}/VEM_Quadrature_3D.cpp)

list(APPEND polydim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_src_sources ${polydim_src_sources} PARENT_SCOPE)
set(polydim_src_headers ${polydim_src_headers} PARENT_SCOPE)
set(polydim_src_includes ${polydim_src_includes} PARENT_SCOPE)
