list(APPEND stokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND stokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)
list(APPEND stokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND stokes_DF_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)

list(APPEND stokes_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND stokes_DF_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND stokes_DF_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(stokes_DF_PCC_2D_sources ${stokes_DF_PCC_2D_sources} PARENT_SCOPE)
set(stokes_DF_PCC_2D_headers ${stokes_DF_PCC_2D_headers} PARENT_SCOPE)
set(stokes_DF_PCC_2D_includes ${stokes_DF_PCC_2D_includes} PARENT_SCOPE)
