list(APPEND navierstokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND navierstokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)
list(APPEND navierstokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND navierstokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)

list(APPEND navierstokes_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND navierstokes_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND navierstokes_DF_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(navierstokes_DF_PCC_2D_sources ${navierstokes_DF_PCC_2D_sources} PARENT_SCOPE)
set(navierstokes_DF_PCC_2D_headers ${navierstokes_DF_PCC_2D_headers} PARENT_SCOPE)
set(navierstokes_DF_PCC_2D_includes ${navierstokes_DF_PCC_2D_includes} PARENT_SCOPE)
