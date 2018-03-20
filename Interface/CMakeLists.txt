# Instructions to compile this program
# are explained in tutorials/README.txt

cmake_minimum_required(VERSION 2.8.11)
project(Interface)

set(CMAKE_PREFIX_PATH "/include/x86_64-linux-gnu/qt5")
#set(Qt5Widgets_DIR /usr/lib/x86_64-linux-gnu/cmake/Qt5Widgets)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
#include(cmake/opengl.cmake)

if (NOT CMAKE_BUILD_TYPE)
  message (STATUS "No build type selected, default to Release")
  set (CMAKE_BUILD_TYPE "Release")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Qt5Core)
find_package(Qt5Widgets)
#find_package(Geogram REQUIRED)

qt5_wrap_cpp(window.cpp renderarea.cpp)

add_executable(interface-build main.cpp renderarea.cpp window.cpp polygoninout.cpp renderarea.h window.h polygoninout.h)

target_link_libraries(interface-build Qt5::Widgets)

target_link_libraries(interface-build voronoi_extractor)

#get_target_property(QtWidgets_location Qt5::Widgets LOCATION)


