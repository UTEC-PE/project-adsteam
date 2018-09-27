cmake_minimum_required(VERSION 3.10)
project(Graph)

set(CMAKE_CXX_STANDARD 17)

add_executable(Graph main.cpp node.h graph.h edge.h read.h)