list(APPEND parabolic_PCC_bulkface_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND parabolic_PCC_bulkface_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND parabolic_PCC_bulkface_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND parabolic_PCC_bulkface_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND parabolic_PCC_bulkface_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND parabolic_PCC_bulkface_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND parabolic_PCC_bulkface_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(parabolic_PCC_bulkface_2D_sources ${parabolic_PCC_bulkface_2D_sources} PARENT_SCOPE)
set(parabolic_PCC_bulkface_2D_headers ${parabolic_PCC_bulkface_2D_headers} PARENT_SCOPE)
set(parabolic_PCC_bulkface_2D_includes ${parabolic_PCC_bulkface_2D_includes} PARENT_SCOPE)
