cmake_minimum_required(VERSION 3.12)

project(Brinkman_DF_PCC_2D LANGUAGES CXX C)

# Create variables
#########################################################

if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()

# Add dependencies
#########################################################

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
      find_package(GeDiM CONFIG REQUIRED)
      list(APPEND Brinkman_DF_PCC_2D_LINKED_LIBRARIES GeDiM::GeDiM)
else()
      get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
      get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
      get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
      get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
      get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

      list(APPEND Brinkman_DF_PCC_2D_LINKED_LIBRARIES ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
      list(APPEND Brinkman_DF_PCC_2D_HEADERS ${GEDIM_LIBRARY_HEADERS})
      list(APPEND Brinkman_DF_PCC_2D_INCLUDE ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## PolyDiM
if (NOT TARGET PolyDiM::PolyDiM)
      find_package(PolyDiM CONFIG REQUIRED)
      list(APPEND Brinkman_DF_PCC_2D_LINKED_LIBRARIES PolyDiM::PolyDiM)
else()
      get_target_property(POLYDIM_LIBRARY_HEADERS PolyDiM::PolyDiM PUBLIC_HEADER)
      get_target_property(POLYDIM_LIBRARY_NAME PolyDiM::PolyDiM EXPORT_NAME)
      get_target_property(POLYDIM_LIBRARY_LINKED_LIBRARIES PolyDiM::PolyDiM LINK_LIBRARIES)
      get_target_property(POLYDIM_LIBRARY_INCLUDE_DIRECTORIES PolyDiM::PolyDiM INCLUDE_DIRECTORIES)
      get_target_property(POLYDIM_LIBRARY_SOURCES PolyDiM::PolyDiM SOURCES)

      list(APPEND Brinkman_DF_PCC_2D_LINKED_LIBRARIES ${POLYDIM_LIBRARY_NAME} ${POLYDIM_LIBRARY_LINKED_LIBRARIES})
      list(APPEND Brinkman_DF_PCC_2D_HEADERS ${POLYDIM_LIBRARY_HEADERS})
      list(APPEND Brinkman_DF_PCC_2D_INCLUDE ${POLYDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## MPI
if(${USE_MPI})
    find_package(MPI REQUIRED)
endif()

## Eigen3
find_package(Eigen3 CONFIG REQUIRED)

## Blas
find_package(BLAS REQUIRED)

## Lapack
find_package(LAPACK REQUIRED)

## VTK
if(${ENABLE_VTK})
      find_package(VTK CONFIG REQUIRED)
endif()

## VORO
if(${ENABLE_VORO})
      find_package(VORO CONFIG REQUIRED)
endif()

## triangle
if(${ENABLE_TRIANGLE})
      find_package(triangle CONFIG REQUIRED)
endif()

## tetgen
if(${ENABLE_TETGEN})
      find_package(tetgen CONFIG REQUIRED)
endif()

## METIS
if(${ENABLE_METIS})
      find_package(metis CONFIG REQUIRED)
endif()

# Insert Sources
#########################################################
add_subdirectory(src)
list(APPEND Brinkman_DF_PCC_2D_HEADERS ${brinkman_DF_PCC_2D_headers})
list(APPEND Brinkman_DF_PCC_2D_SOURCES ${brinkman_DF_PCC_2D_sources} main.cpp)
list(APPEND Brinkman_DF_PCC_2D_INCLUDE ${brinkman_DF_PCC_2D_includes})

# Add variables for clang anc cppcheck
###############################################################################################
list(APPEND CLANGFORMAT_FILES ${brinkman_DF_PCC_2D_sources} "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
list(APPEND CLANGFORMAT_FILES ${brinkman_DF_PCC_2D_headers})
set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)

list(APPEND CPPCHECK_FILES ${brinkman_DF_PCC_2D_sources} "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp")
list(APPEND CPPCHECK_FILES ${brinkman_DF_PCC_2D_headers})
set(CPPCHECK_FILES ${CPPCHECK_FILES} PARENT_SCOPE)

# Create library test executable
#########################################################
add_executable(${PROJECT_NAME}
      ${Brinkman_DF_PCC_2D_SOURCES}
      ${Brinkman_DF_PCC_2D_HEADERS})

target_include_directories(${PROJECT_NAME}
      PRIVATE ${Brinkman_DF_PCC_2D_INCLUDE})

target_link_libraries(${PROJECT_NAME} ${Brinkman_DF_PCC_2D_LINKED_LIBRARIES})

set_target_properties(${PROJECT_NAME} PROPERTIES
      LINKER_LANGUAGE CXX
      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF)
      
# Add test
###############################################################################################
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/integration_test.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
