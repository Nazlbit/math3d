cmake_minimum_required(VERSION 3.0.0)
project(math3d VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

add_library(math3d src/math3d.cpp)

target_include_directories(math3d PUBLIC include)
