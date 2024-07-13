include_directories(../biGraph ../tools)

add_library(pqBicliqeDensest pqBicliqeDensest.h pqBicliqeDensest.cpp rawEdgePivot.h rawEdgePivot.cpp)
target_link_libraries(pqBicliqeDensest biGraph)

add_library(exactFlowAlgorithm exactFlowAlgorithm.h exactFlowAlgorithm.cpp dinic.h rawEdgePivot.h rawEdgePivot.cpp)
target_link_libraries(exactFlowAlgorithm biGraph)
