cmake_minimum_required(VERSION 3.12)

project(PolyDiMUnitTest LANGUAGES CXX C)

# Create library test variables
###############################################################################################

# Add dependencies
###############################################################################################
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

## Threads
find_package(Threads REQUIRED)
list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE Threads::Threads)

## GTest
if (CMAKE_VERSION VERSION_LESS 3.20)
    find_package(GTest REQUIRED)
    find_package(GMock REQUIRED)

    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::GTest)
    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::Main)
    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE ${GMOCK_BOTH_LIBRARIES})
else()
    find_package(GTest REQUIRED)

    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::gtest)
    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::gtest_main)
    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::gmock)
    list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GTest::gmock_main)
endif()

## GeDiM
if (NOT TARGET GeDiM::GeDiM)
	find_package(GeDiM CONFIG REQUIRED)
	list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE GeDiM::GeDiM)
else()
	get_target_property(GEDIM_LIBRARY_HEADERS GeDiM::GeDiM PUBLIC_HEADER)
	get_target_property(GEDIM_LIBRARY_NAME GeDiM::GeDiM EXPORT_NAME)
	get_target_property(GEDIM_LIBRARY_LINKED_LIBRARIES GeDiM::GeDiM LINK_LIBRARIES)
	get_target_property(GEDIM_LIBRARY_INCLUDE_DIRECTORIES GeDiM::GeDiM INCLUDE_DIRECTORIES)
	get_target_property(GEDIM_LIBRARY_SOURCES GeDiM::GeDiM SOURCES)

	list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE ${GEDIM_LIBRARY_NAME} ${GEDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND POLYDIM_LIBRARY_TEST_HEADERS ${GEDIM_LIBRARY_HEADERS})
	list(APPEND POLYDIM_LIBRARY_TEST_INCLUDE ${GEDIM_LIBRARY_INCLUDE_DIRECTORIES})
endif()

## PolyDiM
if (NOT TARGET PolyDiM::PolyDiM)
	find_package(PolyDiM CONFIG REQUIRED)
	list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE PolyDiM::PolyDiM)
else()
	get_target_property(POLYDIM_LIBRARY_HEADERS PolyDiM::PolyDiM PUBLIC_HEADER)
	get_target_property(POLYDIM_LIBRARY_NAME PolyDiM::PolyDiM EXPORT_NAME)
	get_target_property(POLYDIM_LIBRARY_LINKED_LIBRARIES PolyDiM::PolyDiM LINK_LIBRARIES)
	get_target_property(POLYDIM_LIBRARY_INCLUDE_DIRECTORIES PolyDiM::PolyDiM INCLUDE_DIRECTORIES)
	get_target_property(POLYDIM_LIBRARY_SOURCES PolyDiM::PolyDiM SOURCES)

	list(APPEND POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES PRIVATE ${POLYDIM_LIBRARY_NAME} ${POLYDIM_LIBRARY_LINKED_LIBRARIES})
	list(APPEND POLYDIM_LIBRARY_TEST_HEADERS ${POLYDIM_LIBRARY_HEADERS})
	list(APPEND POLYDIM_LIBRARY_TEST_INCLUDE ${POLYDIM_LIBRARY_INCLUDE_DIRECTORIES})
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
list(APPEND POLYDIM_LIBRARY_TEST_PUBLIC_HEADERS ${polydim_test_src_headers})
list(APPEND POLYDIM_LIBRARY_TEST_SOURCES ${polydim_test_src_sources})
list(APPEND POLYDIM_LIBRARY_TEST_SOURCES main.cpp)
list(APPEND POLYDIM_LIBRARY_TEST_HEADERS ${polydim_test_src_headers})
list(APPEND POLYDIM_LIBRARY_TEST_INCLUDE ${polydim_test_src_includes})

# Add variables for clang anc cppcheck
###############################################################################################
list(APPEND CLANGFORMAT_FILES ${POLYDIM_LIBRARY_TEST_SOURCES})
list(APPEND CLANGFORMAT_FILES ${POLYDIM_LIBRARY_TEST_HEADERS})
set(CLANGFORMAT_FILES ${CLANGFORMAT_FILES} PARENT_SCOPE)

list(APPEND CPPCHECK_FILES ${POLYDIM_LIBRARY_TEST_SOURCES})
list(APPEND CPPCHECK_FILES ${POLYDIM_LIBRARY_TEST_HEADERS})
set(CPPCHECK_FILES ${CPPCHECK_FILES} PARENT_SCOPE)

# Create library test executable
###############################################################################################
add_executable(${PROJECT_NAME}
    ${POLYDIM_LIBRARY_TEST_SOURCES}
    ${POLYDIM_LIBRARY_TEST_HEADERS})

target_link_libraries(${PROJECT_NAME} ${POLYDIM_LIBRARY_TEST_LINKED_LIBRARIES})
target_include_directories(${PROJECT_NAME} PRIVATE ${POLYDIM_LIBRARY_TEST_INCLUDE})

set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER "${POLYDIM_LIBRARY_TEST_PUBLIC_HEADERS}")

gtest_discover_tests(${PROJECT_NAME})
