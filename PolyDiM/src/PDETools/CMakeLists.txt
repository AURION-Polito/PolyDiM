add_subdirectory(DOFs)
add_subdirectory(Equations)
add_subdirectory(Assembler)
add_subdirectory(Mesh)
add_subdirectory(LocalSpace)

list(APPEND polydim_src_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(polydim_src_sources ${polydim_src_sources} PARENT_SCOPE)
set(polydim_src_headers ${polydim_src_headers} PARENT_SCOPE)
set(polydim_src_includes ${polydim_src_includes} PARENT_SCOPE)
