list(APPEND SUPG_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND SUPG_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND SUPG_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND SUPG_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND SUPG_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND SUPG_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND SUPG_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(SUPG_PCC_2D_sources ${SUPG_PCC_2D_sources} PARENT_SCOPE)
set(SUPG_PCC_2D_headers ${SUPG_PCC_2D_headers} PARENT_SCOPE)
set(SUPG_PCC_2D_includes ${SUPG_PCC_2D_includes} PARENT_SCOPE)
