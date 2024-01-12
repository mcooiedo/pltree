cmake_minimum_required(VERSION 3.22.1)
project(pltree-ff)

#set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS  "-Iinclude -mclflushopt -mclwb  -march=native -O0 -g  -mavx2 -mrtm ")#-mavx512f -mavx512vl -mavx512bw -mavx512dq -mavx512cd-ltcmalloc
file(GLOB INCLUDE_H "src/include/*.cpp" "src/include/*.h"  "src/include/*.h"  "src/util/*.h"  "src/benchmark/*.cpp" "src/benchmark/*.h" )
include_directories("src/util")
include_directories("src/include")
add_executable(pltree-ff ${INCLUDE_H}
        src/include/fast&fair.cpp)
target_link_libraries(${PROJECT_NAME} pmemobj pthread tbb )  #tcmalloc_minimal
