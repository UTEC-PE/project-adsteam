cmake_minimum_required(VERSION 3.10)
project(Graph)

set(CMAKE_CXX_STANDARD 17)

add_executable(Project main.cpp
                        graph.h
                        disjoint.h
                        )