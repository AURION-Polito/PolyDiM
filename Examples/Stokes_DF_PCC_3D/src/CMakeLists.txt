list(APPEND stokes_DF_PCC_3D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND stokes_DF_PCC_3D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)
list(APPEND stokes_DF_PCC_3D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND stokes_DF_PCC_3D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)

list(APPEND stokes_DF_PCC_3D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND stokes_DF_PCC_3D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND stokes_DF_PCC_3D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(stokes_DF_PCC_3D_sources ${stokes_DF_PCC_3D_sources} PARENT_SCOPE)
set(stokes_DF_PCC_3D_headers ${stokes_DF_PCC_3D_headers} PARENT_SCOPE)
set(stokes_DF_PCC_3D_includes ${stokes_DF_PCC_3D_includes} PARENT_SCOPE)
