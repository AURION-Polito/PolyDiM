list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_2D_LocalSpace.hpp)
list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_2D_Inertia_LocalSpace.hpp)
list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_2D_Ortho_LocalSpace.hpp)

list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_3D_LocalSpace.hpp)
list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_3D_Ortho_LocalSpace.hpp)
list(APPEND polydim_test_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_VEM_PCC_3D_Inertia_LocalSpace.hpp)

list(APPEND polydim_test_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_test_sources ${polydim_test_sources} PARENT_SCOPE)
set(polydim_test_headers ${polydim_test_headers} PARENT_SCOPE)
set(polydim_test_includes ${polydim_test_includes} PARENT_SCOPE)
