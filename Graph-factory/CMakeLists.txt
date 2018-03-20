# Instructions to compile this program
# are explained in tutorials/README.txt

cmake_minimum_required(VERSION 2.8.11)
project(project-library)

include(cmake/opengl.cmake)



if (NOT CMAKE_BUILD_TYPE)
  message (STATUS "No build type selected, default to Release")
  set (CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)


find_package(Geogram REQUIRED)

add_library(voronoi_extractor STATIC Graph.cpp graphmaker.cpp)
#add_executable(geogram_demo_Delaunay2d main.cpp Graph.cpp)

SET_PROPERTY(TARGET voronoi_extractor APPEND PROPERTY
                            INTERFACE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(voronoi_extractor Geogram::geogram Geogram::geogram_gfx  ${GLFW_LIBRARIES})


