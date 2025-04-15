
#include <random>

#include "../biGraph/biGraph.hpp"
#include "../tools/linearSet.hpp"
class accuracy {
   private:
    biGraph* g;
    int p, q;
    LinearSet candL, candR, candRDP;
    std::vector<std::vector<double> > ansAll;
    int minPQ;

    double **C, *bf3;
    void computeC() {
        printf("minPQ: %d\n", minPQ);
        int maxPQ = std::max(minPQ, q - p) + 2;
        int maxD = 5 * std::max(g->maxDu, g->maxDv) + 1;
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

    void computeC2() {
        printf("minPQ: %d\n", minPQ);
        int maxPQ = std::max(minPQ + 1, q - p + 1) + 2;
        int maxD = std::max(g->n1, g->n2) + 1;
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
    ~accuracy() {
        delete g;
        delete[] C;
        delete[] bf3;
    }
    accuracy(const std::string& filePath, const std::string& outFilePath, int p_, int q_, std::string algo = "zstar") {
        p = p_;
        q = q_;
        minPQ = std::min(p, q);
        if (p <= q) {
            g = new biGraph(filePath);
        } else {
            g = new biGraph(filePath, p, q);
        }

        std::printf("load graph\n");
        fflush(stdout);
        if (algo == "zstar4") {
            computeC2();
        } else {
            computeC();
        }
    }

    void testSubgraphSize();
    void approximateCountingAllVersion2(uint32_t T);
    void shadowBuilder1(int p, int q, double e);
    void shadowBuilderZStar(int p, int q, double e);
    void shadowBuilderZStar2(int p, int q, double e);
    void shadowBuilderZStar3(int p, int q, double e,double delta);
    void shadowBuilderZStar4(int p, int q, double e,double delta);
    void shadowBuilderZStar5(int p, int q, double e,double delta);
    void buildDP(int pL, int pR);
    void sampleOne(int length);
    void shadowBuilderAlias(int p, int q, double e);
    int getIndexAlias(std::mt19937& gen, std::vector<double>& Prob, std::vector<uint32_t>& Alias);
    void initializeAliasMethod(std::vector<double>& probabilities, std::vector<double>& Prob, std::vector<uint32_t>& Alias);
    std::vector<uint32_t> reservoirSample(std::vector<uint32_t>& vec, int n);
    std::vector<uint32_t> sampeleRest(int r, int preU, std::vector<uint32_t>& pU, uint32_t outIndex);
    std::vector<uint32_t> sampeleRest2(int r, int preU, std::vector<uint32_t>& pU, uint32_t outIndex);
};
