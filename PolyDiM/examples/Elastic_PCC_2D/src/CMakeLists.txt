list(APPEND elastic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND elastic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND elastic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND elastic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)
list(APPEND elastic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/local_space.hpp)

list(APPEND elastic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND elastic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)
list(APPEND elastic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/local_space.cpp)

list(APPEND elastic_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(elastic_PCC_2D_sources ${elastic_PCC_2D_sources} PARENT_SCOPE)
set(elastic_PCC_2D_headers ${elastic_PCC_2D_headers} PARENT_SCOPE)
set(elastic_PCC_2D_includes ${elastic_PCC_2D_includes} PARENT_SCOPE)
