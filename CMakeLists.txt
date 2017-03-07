cmake_minimum_required(VERSION 3.6)
project(migra_gapbs)

set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_COMPILER /Users/yunmingzhang/Documents/intel/bin/icpc)
set(CMAKE_CXX_FLAGS "-openmp")

include_directories("./src")

add_executable(bc src/bc.cc)
add_executable(bfs src/bfs.cc)
add_executable(me src/multi_hop_simple_rec.cc)