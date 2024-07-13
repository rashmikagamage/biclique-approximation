# This is the source code of "Efficient Biclique Counting in Large Bipartite Graphs, SIGMOD_2023" by Xiaowei Ye

EPivoter:
./run -f datafile -pm

Zigzag
./run -f datafile -cp -t 100000

Zigzag++
./run -f datafile -cp -t 100000 -v5

EP/Zigzag
./run -f datafile -pp -t 100000

EP/Zigzag++
./run -f datafile -pp -t 100000 -v5

Other usage details is in the run.cpp.

VLDB17:the implementation of "On Sampling from Massive Graph Streams" for biclique counting

src/densestSubgraph/: the implementation of the (p,q)-biclique densest subgraph problem

to get the higher-order clustering coefficient, uncomment the line 13 "#define PQWEDGE" of src/biCliique/rawEdgePivot.h and recompile. Note the the output ranges from 0 to 0.5. You need to double the output.
# biclique-approximation
# biclique-approximation
