list(APPEND brinkman_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND brinkman_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)
list(APPEND brinkman_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND brinkman_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)

list(APPEND brinkman_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND brinkman_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND brinkman_DF_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(brinkman_DF_PCC_2D_sources ${brinkman_DF_PCC_2D_sources} PARENT_SCOPE)
set(brinkman_DF_PCC_2D_headers ${brinkman_DF_PCC_2D_headers} PARENT_SCOPE)
set(brinkman_DF_PCC_2D_includes ${brinkman_DF_PCC_2D_includes} PARENT_SCOPE)
