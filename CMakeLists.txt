cmake_minimum_required(VERSION 2.8.11)
project(Voronoi-All)

if (NOT CMAKE_BUILD_TYPE)
  message (STATUS "No build type selected, default to Release")
  set (CMAKE_BUILD_TYPE "Release")
#  set (CMAKE_BUILD_TYPE "Debug")
endif ()

add_subdirectory(Clipper)
add_subdirectory(Graph-factory)
add_subdirectory(MinimalTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
