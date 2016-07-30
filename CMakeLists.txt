cmake_minimum_required(VERSION 3.5)
project(SPH)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES main.cpp Vector2.h)
add_executable(SPH ${SOURCE_FILES})