cmake_minimum_required(VERSION 2.8)
project(RapidMatch)

include_directories(
        ${PROJECT_SOURCE_DIR}/utility
        ${PROJECT_SOURCE_DIR}/
)


set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -std=c++14 -O3 -g -Wall -pthread -march=native")

add_subdirectory(graph)
add_subdirectory(matching)
add_subdirectory(utility)
add_subdirectory(preprocessing)
