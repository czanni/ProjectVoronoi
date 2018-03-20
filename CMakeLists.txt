cmake_minimum_required(VERSION 2.8.11)
project(Voronoi-All)

add_subdirectory(Graph-factory)
add_subdirectory(Interface)


set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
