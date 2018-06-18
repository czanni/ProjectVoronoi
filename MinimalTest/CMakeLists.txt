# Instructions to compile this program
# are explained in tutorials/README.txt

cmake_minimum_required(VERSION 2.8.11)
project(MinimalTest)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)

add_executable(graphic-build main.cpp polygoninout.cpp)

target_link_libraries(graphic-build voronoi_extractor)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
#include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )
target_link_libraries(graphic-build ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )

