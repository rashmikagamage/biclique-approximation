
#include <random>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"

struct Node {
    std::vector<uint32_t> SU;
    std::vector<uint32_t> SV;
    std::vector<Node> children;
    int label;

    Node(const std::vector<uint32_t>& su = {},
         const std::vector<uint32_t>& sv = {},
         int lbl = 0)
        : SU(su), SV(sv), label(lbl) {}
};

class BCT {
   private:
    biGraph* g;
    int p, q;
    LinearSet candL, candR;
    std::vector<std::vector<double> > ansAll;
    int minPQ;

    double **C, *bf3;
    void computeC() {
        printf("minPQ: %d\n", minPQ);
        int maxPQ = std::max(minPQ, q - p) + 2;
        int maxD = std::max(g->maxDu, g->maxDv) + 1;
        C = new double*[maxD];
        bf3 = new double[maxD * maxPQ];
        for (int i = 0; i < maxD; i++) {
            C[i] = bf3 + i * maxPQ;
        }
        C[0][0] = 1;
        C[1][0] = 1;
        C[1][1] = 1;
        for (int i = 2; i < maxD; i++) {
            C[i][0] = 1;
            if (i < maxPQ) C[i][i] = 1;
            for (int j = 1; j < i && j < maxPQ; j++) {
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
    std::pair<uint32_t, bool> selectPivot(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV);

    biGraph createSubgraph(const std::vector<uint32_t>& SU, const std::vector<uint32_t>& SV);
    void printPath(const std::vector<int>& path);
    void dfs(Node& node, std::vector<int>& currentPath);
    void traversePaths(Node& root);
};