#ifndef COLORPATH
#define COLORPATH

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

class colorPath {
    
private:

    biGraph * g;
    int p, q;
    LinearSet candL, candR;
    std::vector<std::vector<double> > ansAll;
    int minPQ;

    double ** C, *bf3;
    void computeC() {
        
        // int maxPQ = std::max(p, q) + 1;
        int maxPQ = minPQ + 1;
        int maxD = std::max(g->maxDu, g->maxDv) + 1;

        C = new double*[maxD];
        bf3 = new double[maxD * maxPQ];
        for(int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxPQ;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for(int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if(i < maxPQ) C[i][i] = 1;
            for(int j = 1; j < i && j < maxPQ; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }

public:
    ~colorPath() {
        delete g;
        delete [] C;
        delete [] bf3;
    }
    colorPath(const std::string & filePath, const std::string & outFilePath, 
        int H=1, int p_=0, int q_=0) {
        p = p_;
        q = q_;
        minPQ = H;

        g = new biGraph(filePath);
        printf("load graph\n");fflush(stdout);

        computeC();
    }

    void approximateCounting(uint64_t T);
    void testSubgraphSize();
    void approximateCounting2(uint64_t T);
    void approximateCounting3(uint64_t T);
    
    void approximateCounting4(uint64_t T);

    void approximateCountingAll(uint64_t T);
//差分优化
    void approximateCountingAllVersion2(uint64_t T);
//点层面分割
    void approximateCountingAllVersion3(uint64_t T);
//整个图dp
    void approximateCountingAllVersion4(uint64_t T);
//快速建图
    void approximateCountingAllVersion5(uint64_t T);

    void buildDP(int minPQ);

    void buildShadowFixed(int p,int q);
    void buildShadow1(int p, int  q , double epsilon);
};

#endif