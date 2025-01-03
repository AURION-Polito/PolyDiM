cmake_minimum_required(VERSION 3.12)
project(PolyDiM_Examples LANGUAGES CXX C VERSION 1.0)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_PCC_1D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_PCC_2D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_PCC_3D)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_MCC_2D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_MCC_3D)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Brinkman_DF_PCC_2D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/NavierStokes_DF_PCC_2D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Stokes_DF_PCC_3D)

set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)
