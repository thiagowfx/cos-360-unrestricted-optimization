cmake_minimum_required(VERSION 2.8.12)
set(PROJECT_NAME "main")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O2 -std=c++11")
project(${PROJECT_NAME} C CXX)

set(EXT_PROJECTS_DIR "${PROJECT_SOURCE_DIR}/ext")
set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
set(TEST_DIR "${PROJECT_SOURCE_DIR}/test")

set(
  PROJECT_SOURCES
  "${SRC_DIR}/lib.hpp"
  )

include_directories(
  ${SRC_DIR}
  )

add_library(
  PROJECT_LIBRARY 
  OBJECT
  ${PROJECT_SOURCES}
  )

add_executable(
  ${PROJECT_NAME}
  ${SRC_DIR}/main.cpp
  $<TARGET_OBJECTS:PROJECT_LIBRARY>
  )

option(TEST "Build all tests." ON)
if (TEST)
  enable_testing()
  set(PROJECT_TEST_NAME "tests")
  file(GLOB TEST_SRC_FILES ${TEST_DIR}/*Test.cpp)
  add_subdirectory(${EXT_PROJECTS_DIR}/gtest)
  include_directories(${GTEST_INCLUDE_DIRS})
  add_executable(
    ${PROJECT_TEST_NAME}
    ${TEST_SRC_FILES}
    $<TARGET_OBJECTS:PROJECT_LIBRARY>
    )
  add_dependencies(${PROJECT_TEST_NAME} googletest)
  target_link_libraries(
    ${PROJECT_TEST_NAME}
    ${GTEST_LIBS_DIR}/libgtest.a
    ${GTEST_LIBS_DIR}/libgtest_main.a
    pthread
    )
  GTEST_ADD_TESTS(${PROJECT_TEST_NAME} "" ${TEST_SRC_FILES})
endif()
