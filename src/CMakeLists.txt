cmake_minimum_required(VERSION 3.6)

project(run)

SET(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -march=native -mavx -O3")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(exName1 "run")
set(exName2 "densest")
set(exName3 "exdensest")

add_subdirectory(biClique)
add_subdirectory(biGraph)
add_subdirectory(tools)
add_subdirectory(densestSubgraph)

add_executable(${exName1} runner/run.cpp)
add_executable(${exName2} runner/densestSubgraph.cpp)
add_executable(${exName3} runner/exactDensestSubgraph.cpp)

target_link_libraries(${exName1} getArgs)
target_link_libraries(${exName1} biGraph)
target_link_libraries(${exName1} BCT)
target_link_libraries(${exName1} BCListPlusPlus)
target_link_libraries(${exName1} BK)
target_link_libraries(${exName1} rawEdgePivot)
target_link_libraries(${exName1} fastEdgePivot)
target_link_libraries(${exName1} edgePivotSpecificPQ)
target_link_libraries(${exName1} colorPath)
target_link_libraries(${exName1} accuracy)

target_link_libraries(${exName1} pivotAndPath)
target_link_libraries(${exName1} turan)
target_link_libraries(${exName1} bcAndPath)
target_link_libraries(${exName1} colorPathPequalsQ)
target_link_libraries(${exName1} pivotAndPathPequalsQ)
target_link_libraries(${exName1} colorPathSpecificPQ)

target_link_libraries(${exName2} getArgs)
target_link_libraries(${exName2} biGraph)
target_link_libraries(${exName2} pqBicliqeDensest)

target_link_libraries(${exName3} getArgs)
target_link_libraries(${exName3} biGraph)
target_link_libraries(${exName3} exactFlowAlgorithm)


# set(binPath "./bin")
# ADD_CUSTOM_COMMAND(TARGET ${exName1} ${exName2} ${exName3}
#           POST_BUILD
#           COMMAND ${CMAKE_COMMAND} -E make_directory ${binPath}
#           COMMAND ${CMAKE_COMMAND} -E copy ${exName1} ${binPath}/
#           COMMAND ${CMAKE_COMMAND} -E copy ${exName2} ${binPath}/
# )
