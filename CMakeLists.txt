cmake_minimum_required(VERSION 3.5)
project(circular_stat)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    circular_stat.cpp
    circular_stat.h
    main.cpp)

add_executable(circular_stat ${SOURCE_FILES})