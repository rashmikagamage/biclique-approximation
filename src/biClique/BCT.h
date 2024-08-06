
#include <algorithm>
#include <random>
#include <utility>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
#include "unordered_set"
#include "vector"

struct Node {
    std::unordered_set<uint32_t> SU;
    std::unordered_set<uint32_t> SV;
    std::vector<Node*> children;
    int vertexSide;  // u is 1 v is 2
    int label;
    int vertex;
    int pCount;
    int vh = 0, vp = 0, uh = 0, up = 0;

    Node(const std::unordered_set<uint32_t>& su = {},
         const std::unordered_set<uint32_t>& sv = {},
         int lbl = 0,
         int vertexSide = 0,
         int vertex = 0,
         int pC = 0)  // u is 1 v is 2
        : SU(su), SV(sv), label(lbl), vertexSide(vertexSide), vertex(vertex), pCount(pC) {}

    ~Node() {
        for (Node* child : children) {
            delete child;
        }
    }
};

class BCT {
   private:
    biGraph* g;
    int p, q;
    LinearSet candL, candR;
    std::vector<std::vector<double>> ansAll;
    int minPQ;
    std::vector<std::vector<double>> count;
    double **C, *bf3;
    void computeC() {
        int maxPQ = g->maxDu * g->maxDv + 1;
        int maxC = std::min(maxPQ, 2000);
        C = new double*[maxC];
        bf3 = new double[maxC * maxC];

        for (int i = 0; i < maxC; i++) {
            C[i] = bf3 + i * maxC;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;

        for (int i = 2; i < maxC; i++) {
            C[i][0] = 1;
            if (i < maxC) C[i][i] = 1;
            for (int j = 1; j < i && j < maxC; j++) {
                C[i][j] = C[i - 1][j - 1] + C[i - 1][j];
            }
        }
    }

   public:
    ~BCT() {
        delete g;
        delete[] C;
        delete[] bf3;
    }
    BCT(const std::string& filePath, const std::string& outFilePath, int p_, int q_) {
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
    std::pair<uint32_t, bool> selectPivot(std::unordered_set<uint32_t>& SU, std::unordered_set<uint32_t>& SV);
    std::pair<uint32_t, bool> selectPivoteWithSide(std::unordered_set<uint32_t>& SU, std::unordered_set<uint32_t>& SV, int pivotSide);
    biGraph createSubgraph(const std::vector<uint32_t>& SU, const std::vector<uint32_t>& SV);
    void printPath(const std::vector<int>& path);
    void dfs(Node& node, std::vector<int>& currentPath);
    void traversePaths(Node* node, std::vector<Node*>& path, std::vector<std::vector<Node*>>& allPaths);
    std::vector<std::vector<double>> countBicliques(Node* root);
    double countBicliquesForPQ(Node* root);
    void processNode(Node* node);
};