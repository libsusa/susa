cmake_minimum_required(VERSION 2.8)
project(susa)

set(BUILD_SHARED_LIBS ON)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

if(UNIX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11")
endif(UNIX)

if (${OPTIMIZED_BUILD})
  message(STATUS "Optimized build is ON")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(${OPTIMIZED_BUILD})

include_directories(inc)
set(SRC_FILES src/rng.cpp
              src/mt.cpp
              src/base.cpp
              src/convolutional.cpp
              src/matrix.cpp
              src/modulation.cpp
              src/rrcosine.cpp
              src/svd.cpp
              src/utility.cpp)

add_library(susa ${SRC_FILES})

add_subdirectory (test)
enable_testing ()
add_test (NAME check COMMAND ./test/check)
