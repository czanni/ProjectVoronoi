cmake_minimum_required(VERSION 2.8.11)
project(project-library)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

add_library(voronoi_extractor STATIC Graph.cpp GraphMaker.cpp Delaunay_psm.cpp)

SET_PROPERTY(TARGET voronoi_extractor APPEND PROPERTY
                            INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(voronoi_extractor Threads::Threads)

target_link_libraries(voronoi_extractor polyclipping)

# To compile Delaunay_psm under Linux:
#   g++ -O3 -fopenmp -frounding-math -ffp-contract=off Delaunay_example.cpp Delaunay_psm.cpp -o Delaunay_example -ldl -lm
#TODO a few problem with the original file
#  templated function in psm.cpp
#  different signature for triangle_circumcenter
#  require "using isnan;" for other plateform than mac
