cmake_minimum_required(VERSION 3.6)

project(run)

# Use g++ as the compiler
SET(CMAKE_CXX_COMPILER "g++")

# Set C++11 standard and optimization flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mavx -O0 -g")

# Ensure asserts are always enabled by undefining NDEBUG (even if defined in code)
add_definitions(-U NDEBUG)

# Set output directory for binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Define executable names
set(exName1 "run")
set(exName2 "densest")
set(exName3 "exdensest")

# Include subdirectories
add_subdirectory(biClique)
add_subdirectory(biGraph)
add_subdirectory(tools)
add_subdirectory(densestSubgraph)

# Add executables
add_executable(${exName1} runner/run.cpp)
add_executable(${exName2} runner/densestSubgraph.cpp)
add_executable(${exName3} runner/exactDensestSubgraph.cpp)

# Link necessary libraries for each executable
target_link_libraries(${exName1} getArgs biGraph BCListPlusPlus BK rawEdgePivot fastEdgePivot edgePivotSpecificPQ colorPath accuracy pivotAndPath turan bcAndPath colorPathPequalsQ pivotAndPathPequalsQ colorPathSpecificPQ)
target_link_libraries(${exName2} getArgs biGraph pqBicliqeDensest)
target_link_libraries(${exName3} getArgs biGraph exactFlowAlgorithm)
