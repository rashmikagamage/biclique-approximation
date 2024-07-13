cmake_minimum_required(VERSION 3.6)

project(run)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mavx -O3")

set(exName1 "run")

add_executable(${exName1} main.cpp)

add_library(fastIO INTERFACE)
add_library(getArgs  INTERFACE)
add_library(hopstotchHash INTERFACE)
add_library(linearSet INTERFACE)
add_library(listLinearHeap  INTERFACE)
add_library(biGraph INTERFACE)

add_library(BCListPlusPlus BCListPlusPlus.h BCListPlusPlus.cpp)
target_link_libraries(BCListPlusPlus biGraph linearSet)

target_link_libraries(${exName1} getArgs)
target_link_libraries(${exName1} biGraph)
target_link_libraries(${exName1} BCListPlusPlus)

