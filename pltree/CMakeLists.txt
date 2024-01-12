cmake_minimum_required(VERSION 3.22.1)
project(pltree)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS  "-Iinclude -mclflushopt -mclwb  -march=native -O3 -mavx512f -mavx512vl -mavx512bw -mavx512dq -mavx512cd -g  -mavx2 ")#-ltcmalloc
file(GLOB INCLUDE_H "src/include/*.cpp" "src/include/*.h"  "src/include/*.h"  "src/util/*.h"  "src/benchmark/*.cpp" "src/benchmark/*.h" )
include_directories("src/util")
include_directories("src/include")
add_executable(pltree ${INCLUDE_H} )
target_link_libraries(${PROJECT_NAME} pmemobj pthread tbb )  #tcmalloc_minimal
