cmake_minimum_required(VERSION 3.2)
project(susa)
enable_testing()


if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")


include_directories (inc)
set (SRC_FILES src/rng.cpp
               src/mt.cpp
               src/sets.cpp
               src/base.cpp
               src/convolutional.cpp
               src/matrix.cpp
               src/modulation.cpp
               src/rrcosine.cpp
               src/svd.cpp
               src/utility.cpp)

add_library (susa SHARED ${SRC_FILES})

add_subdirectory (examples)

add_subdirectory (test)

install (TARGETS susa
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib)

install (DIRECTORY inc/
         DESTINATION include
         FILES_MATCHING PATTERN "*.h")
