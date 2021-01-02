cmake_minimum_required(VERSION 3.2)
project(susa)
enable_testing()


if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif(APPLE)

set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror -pedantic")

option(CODE_COVERAGE "Enable coverage reporting" OFF)
if(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  # Add required flags (GCC & LLVM/Clang)
  add_compile_options(
    -O0        # no optimization
    -g         # generate debug info
    --coverage # sets all required flags
  )

  link_libraries(--coverage)
endif(CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")


include_directories (inc)
set (SRC_FILES src/rng.cpp
               src/mt.cpp
               src/sets.cpp
               src/base.cpp
               src/ccode.cpp
               src/matrix.cpp
               src/modulation.cpp
               src/rrcosine.cpp
               src/svd.cpp
               src/utility.cpp)

add_library (susa SHARED ${SRC_FILES})

option(EXAMPLES "Enable examples build" ON)
if (EXAMPLES)
add_subdirectory (examples)
endif(EXAMPLES)

add_subdirectory (test)

install (TARGETS susa
         LIBRARY DESTINATION lib
         ARCHIVE DESTINATION lib)

install (DIRECTORY inc/
         DESTINATION include
         FILES_MATCHING PATTERN "*.h")
