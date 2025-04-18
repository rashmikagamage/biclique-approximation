cmake_minimum_required(VERSION 3.16)

# Set the C++ compiler and standard
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set compiler flags for optimizations and hardware features
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mavx -mavx2 -O0 -g")

# Always enable assertions by undefining NDEBUG in all build types
add_definitions(-U NDEBUG)

# Include necessary directories
include_directories(../biGraph ../tools)

# Add libraries and link them as needed



add_library(BCListPlusPlus BCListPlusPlus.h BCListPlusPlus.cpp)
target_link_libraries(BCListPlusPlus biGraph linearSet)

add_library(BK BK.h BK.cpp)
target_link_libraries(BK biGraph linearSet)

add_library(rawEdgePivot rawEdgePivot.h rawEdgePivot.cpp)
target_link_libraries(rawEdgePivot biGraph linearSet)

add_library(fastEdgePivot fastEdgePivot.h fastEdgePivot.cpp)
target_link_libraries(fastEdgePivot biGraph linearSet)

add_library(edgePivotSpecificPQ edgePivotSpecificPQ.h edgePivotSpecificPQ.cpp)
target_link_libraries(edgePivotSpecificPQ biGraph linearSet)

add_library(colorPath colorPath.h colorPath.cpp colorPathAV5.cpp)
target_link_libraries(colorPath biGraph linearSet)

add_library(pivotAndPath pivotAndPath.h pivotAndPath.cpp pivotAndPathSamplev5.cpp)
target_link_libraries(pivotAndPath biGraph linearSet)

add_library(turan turan.h turan.cpp)
target_link_libraries(turan biGraph)

add_library(bcAndPath bcAndPath.h bcAndPath.cpp bcAndPathV5.cpp)
target_link_libraries(bcAndPath biGraph linearSet)

add_library(colorPathPequalsQ colorPathPequalsQ.h colorPathPequalsQ.cpp colorPathAV5PequalsQ.cpp)
target_link_libraries(colorPathPequalsQ biGraph linearSet)

add_library(pivotAndPathPequalsQ pivotAndPathPequalsQ.h pivotAndPathPequalsQ.cpp pivotAndPathSamplev5PequalsQ.cpp)
target_link_libraries(pivotAndPathPequalsQ biGraph linearSet)

add_library(colorPathSpecificPQ colorPathSpecificPQ.h colorPathSpecificPQ.cpp)
target_link_libraries(colorPathSpecificPQ biGraph linearSet)

add_library(accuracy accuracy.h accuracy.cpp)
target_link_libraries(accuracy biGraph linearSet)

