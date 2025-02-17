cmake_minimum_required(VERSION 3.12)

project(PolyDiM LANGUAGES CXX C VERSION 1.0.0)

# Create library variables
###################################################################
set(POLYDIM_LIBRARY_COMPILE_DEFINITIONS PolyDiM_compile_definitions)
set(POLYDIM_INSTALL_DIR ${PROJECT_NAME})
set(POLYDIM_INSTALL_BINARY_DIR bin/)
set(POLYDIM_INSTALL_INCLUDE_DIR include/)
set(POLYDIM_INSTALL_LIB_DIR lib/)
set(POLYDIM_INSTALL_ARCHIVE_DIR lib/)
set(POLYDIM_INSTALL_DESTINATION lib/cmake/${POLYDIM_INSTALL_DIR})
###################################################################

if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel" FORCE)
endif()
Message(STATUS "PolyDiM Build configuration: " ${CMAKE_BUILD_TYPE})

# IMPOSE WARNINGS ON DEBUG
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall -Wextra -pedantic-errors")

# IMPOSE CXX FLAGS FOR WINDOWS
if (WIN32)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-mbig-obj")
endif (WIN32)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/${POLYDIM_INSTALL_DIR} CACHE PATH "Install path prefix, prepended onto install directories." FORCE)
Message(STATUS "PolyDiM Library will be installed in: " ${CMAKE_INSTALL_PREFIX})


# Add dependencies
###########################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
      find_package(GeDiM CONFIG REQUIRED)
      list(APPEND POLYDIM_LIBRARY_LINKED_LIBRARIES PRIVATE GeDiM::GeDiM)
else()
      get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
      get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
      get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
      get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
      get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

      list(APPEND POLYDIM_LIBRARY_LINKED_LIBRARIES PRIVATE ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
      list(APPEND POLYDIM_LIBRARY_HEADERS ${GEDIM_LIBRARY_HEADERS})
      list(APPEND POLYDIM_LIBRARY_INCLUDE ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
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
###############################################################################################
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/src)
list(APPEND POLYDIM_LIBRARY_SOURCES ${polydim_src_sources})
list(APPEND POLYDIM_LIBRARY_HEADERS ${polydim_src_headers})
list(APPEND POLYDIM_LIBRARY_INCLUDE ${polydim_src_includes})

# Add variables for clang anc cppcheck
###############################################################################################
list(APPEND CLANGFORMAT_FILES ${polydim_src_sources})
list(APPEND CLANGFORMAT_FILES ${polydim_src_headers})
set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)

list(APPEND CPPCHECK_FILES ${polydim_src_sources})
list(APPEND CPPCHECK_FILES ${polydim_src_headers})
set(CPPCHECK_FILES ${CPPCHECK_FILES} PARENT_SCOPE)

# Create public headers
###############################################################################################
set(POLYDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH ${CMAKE_CURRENT_BINARY_DIR}/PolyDiM_additional_include)

foreach(header ${POLYDIM_LIBRARY_HEADERS})
      get_filename_component(headerDirectory ${header} DIRECTORY)
      get_filename_component(headerFileName ${header} NAME_WE)
      get_filename_component(headerExtension ${header} EXT)

      if(${headerExtension} MATCHES ".in")
            set(headerExtension ".hpp")
            set(headerDirectory ${POLYDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

            configure_file(${header} ${headerDirectory}/${headerFileName}${headerExtension} @ONLY)
      elseif(NOT ${headerExtension} MATCHES ".hpp" AND NOT ${headerExtension} MATCHES ".h")
            message(SEND_ERROR "Header extension not recognized for file ${headerFileName}${headerExtension}")
      endif()

      list(APPEND POLYDIM_LIBRARY_PUBLIC_HEADERS ${headerDirectory}/${headerFileName}${headerExtension})
endforeach()

list(APPEND POLYDIM_LIBRARY_INCLUDE ${POLYDIM_LIBRARY_ADDITIONAL_INCLUDE_PATH})

# Create library
###############################################################################################
add_library(${PROJECT_NAME} STATIC ${POLYDIM_LIBRARY_SOURCES} ${POLYDIM_LIBRARY_HEADERS})
add_library(PolyDiM::PolyDiM ALIAS ${PROJECT_NAME})

target_link_libraries(${PROJECT_NAME} ${POLYDIM_LIBRARY_LINKED_LIBRARIES})

target_include_directories(${PROJECT_NAME} PRIVATE ${POLYDIM_LIBRARY_INCLUDE})
target_include_directories(${PROJECT_NAME} PUBLIC $<INSTALL_INTERFACE:${POLYDIM_INSTALL_INCLUDE_DIR}>)
target_compile_definitions(${PROJECT_NAME} PUBLIC ${${POLYDIM_LIBRARY_COMPILE_DEFINITIONS}})

target_compile_options(${PROJECT_NAME} PUBLIC -fPIC)

# Create install target
###############################################################################################
set_target_properties(${PROJECT_NAME} PROPERTIES 
      PUBLIC_HEADER "${POLYDIM_LIBRARY_PUBLIC_HEADERS}"
      EXPORT_NAME PolyDiM
      LINKER_LANGUAGE CXX
      CXX_STANDARD 20
      CXX_STANDARD_REQUIRED ON
      CXX_EXTENSIONS OFF)

install(
      TARGETS ${PROJECT_NAME}
      EXPORT PolyDiMTargets
      RUNTIME DESTINATION ${POLYDIM_INSTALL_BINARY_DIR} # runtime executable
      LIBRARY DESTINATION ${POLYDIM_INSTALL_LIB_DIR} # dynamic libraries
      ARCHIVE DESTINATION ${POLYDIM_INSTALL_ARCHIVE_DIR} # static libraries
      PUBLIC_HEADER DESTINATION ${POLYDIM_INSTALL_INCLUDE_DIR} # headers
)


include(CMakePackageConfigHelpers)

configure_package_config_file(
      ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PolyDiMConfig.cmake.in
      ${CMAKE_CURRENT_BINARY_DIR}/PolyDiMConfig.cmake
      INSTALL_DESTINATION ${CMAKE_INSTALL_PREFIX}/${POLYDIM_INSTALL_DESTINATION}
      PATH_VARS POLYDIM_INSTALL_INCLUDE_DIR
)

write_basic_package_version_file(
      ${CMAKE_CURRENT_BINARY_DIR}/PolyDiMConfigVersion.cmake
      VERSION ${PolyDiM_VERSION}
      COMPATIBILITY AnyNewerVersion
)

export(
      TARGETS ${PROJECT_NAME}
      NAMESPACE PolyDiM::
      FILE ${CMAKE_CURRENT_BINARY_DIR}/PolyDiMTargets.cmake
)

install(
      EXPORT PolyDiMTargets
      NAMESPACE PolyDiM::
      DESTINATION ${CMAKE_INSTALL_PREFIX}/${POLYDIM_INSTALL_DESTINATION}
)

install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/PolyDiMConfig.cmake
      ${CMAKE_CURRENT_BINARY_DIR}/PolyDiMConfigVersion.cmake
      DESTINATION
      ${CMAKE_INSTALL_PREFIX}/${POLYDIM_INSTALL_DESTINATION}
)
