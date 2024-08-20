
#include <algorithm>
#include <random>
#include <utility>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
#include "vector"

struct NodeV {
    int vertexSide;  // u is 1 v is 2
    int label;
    int vertex;
    int pCount;
    int vh = 0, vp = 0, uh = 0, up = 0;
    int depth = 0;

    NodeV(
        int lbl = 0,
        int vertexSide = 0,
        int vertex = 0,
        int pC = 0)  // u is 1 v is 2
        : label(lbl), vertexSide(vertexSide), vertex(vertex), pCount(pC) {}
};

class BCTV2 {
   private:
    biGraph* g;
    int p, q;
    LinearSet candL, candR;
    int minPQ;
    std::vector<std::vector<double>> count;
    double **C, *bf3;

    ui timePivot = 0;
    int uh = 0, vh = 0, up = 0, vp = 0;
    void computeC() {
        int maxI = std::max(g->maxDu, g->maxDv);
        int maxJ = maxI;

        C = new double*[maxI + 1];
        for (int i = 0; i <= maxI; i++) {
            int rowSize = std::min(i, maxJ) + 1;
            C[i] = new double[rowSize];
        }

        C[0][0] = 1;
        for (int i = 1; i <= maxI; i++) {
            C[i][0] = 1;
            if (i <= maxJ) {
                C[i][i] = 1;
            }
        }
        for (int i = 2; i <= maxI; i++) {
            for (int j = 1; j <= std::min(i - 1, maxJ); j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }

   public:
    ~BCTV2() {
        delete g;
        delete[] C;
        delete[] bf3;
    }
    BCTV2(const std::string& filePath, const std::string& outFilePath, int p_, int q_) {
        p = p_;
        q = q_;
        minPQ = std::min(p, q);
        g = new biGraph(filePath);

        std::printf("load graph\n");
        fflush(stdout);
        computeC();
    }

    void buildTree();
    void buildTreeV2(int p, int q);
    std::pair<uint32_t, bool> selectPivot(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV);
    std::pair<uint32_t, bool> selectPivoteWithSide(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV, int pivotSide);
    biGraph createSubgraph(const std::vector<uint32_t>& SU, const std::vector<uint32_t>& SV);
    void printPath(const std::vector<int>& path);
    void dfs(NodeV& node, std::vector<int>& currentPath);
    void traversePaths(NodeV* node, std::vector<NodeV*>& path, std::vector<std::vector<NodeV*>>& allPaths);
    std::vector<std::vector<double>> countBicliques(NodeV* root);
    double countBicliquesForPQ(NodeV* root);
    void processNode(NodeV* node);
    void destroyVector(std::vector<int>& vec);
};