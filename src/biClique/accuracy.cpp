#include "accuracy.h"
#undef NDEBUG
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <csignal>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../biGraph/biGraph.hpp"
using namespace std;
static uint32_t seg_value = 0;
static int pU_seg = -1;
static int pU1_seg = -1;
// Signal handler for segmentation fault
void segfault_handler(int signal) {
    printf("Segmentation fault. Signal: %d\n", signal);
    printf("seg_value: %d\n", seg_value);
    printf("pU_seg: %d\n", pU_seg);
    printf("pU1_seg: %d\n", pU1_seg);
    std::exit(signal);  // Exit the program after printing
}
struct Subspace {
    int p_prime, q_prime;
    std::vector<uint32_t> SU, SV;
    uint32_t SUSize, SVSize;
    double mu;
    double ZCount;
    int subgraph = 0;
    Subspace(int p_prime, int q_prime, std::vector<uint32_t> SU, std::vector<uint32_t> SV,
             uint32_t SUSize, uint32_t SVSize, double mu, double countPQMean,
             double ZCount)
        : p_prime(p_prime), q_prime(q_prime), SU(SU), SV(SV), SUSize(SUSize), SVSize(SVSize), mu(mu), ZCount(ZCount) {}
};

struct CompareSubspace {
    bool operator()(const Subspace& a, const Subspace& b) { return a.mu > b.mu; }
};

void accuracy::shadowBuilderZStar(int p, int q, double e) {
    if (p > q) {
        swap(p, q);
    }

    printf(" Working on p q e: %d - %d - %.2f \n", p, q, e);
    auto start = chrono::high_resolution_clock::now();
    minPQ = min(p, q);
    // get the total structure count in the input graph
    // e1 hols u's neighbors

    std::vector<uint32_t> candLtemp(g->maxDv);
    std::vector<uint32_t> outPosition(g->maxDu);
    vector<int> twoHopCount(g->n1);
    candL.resize(g->n1);
    candR.resize(g->n2);
    double totalZinShadow = 0;
    double totalZinShadowTest = 0;

    uint32_t maxE = std::min((1u << 21), g->maxDu * g->maxDv);
    // for Shadow
    uint32_t cnt_PQ = 1;

    priority_queue<Subspace, std::vector<Subspace>, CompareSubspace> shadow;
    vector<Subspace> subspaceTest;
    uint32_t totalS = 0;
    // DP Building

    std::vector<uint32_t> pV(g->maxDu);  // Adjusted size to accommodate pV[v + 1]
    std::vector<uint32_t> pU(g->maxDv);
    std::vector<uint32_t> e1(maxE);  // Assuming maximum possible size
    std::vector<uint32_t> e2(
        maxE);  // Assuming maximum possible sizmaxPossibleLengthe
    std::vector<double> ddp(maxE, 0.0);
    std::vector<std::vector<double>> dpU(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<std::vector<double>> dpV(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<uint32_t> mapUtoV(maxE);
    std::vector<uint32_t> mapVtoU(maxE);
    const double epsilon = 0.05;
    const double delta = 0.05;
    double gamma = 4 * (1 + epsilon) * (std::exp(1) - 2) * std::log(2 / delta) /
                   (epsilon * epsilon);

    double tSample = 0;
    double tTotal = 0;
    double auxTotalSamples = 1000;
    double pqBicluqeEstAux = 1;
    int candLSum = 0;
    int candRSum = 0;

    for (uint32_t e = 0; e < maxE; e++) {
        dpU[0][e] = 0;
        dpU[1][e] = 1;
        dpV[0][e] = 1;
    }

    // DP Buidling
    auto buildDP = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            if (rSize < g->deg1(x)) {
                for (int k = 0; k < rSize; k++) {
                    uint32_t y = candR[k];
                    if (g->connectUV(x, y)) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            } else {
                auto start = g->e1.begin() + g->pU[x];
                auto end = g->e1.begin() + g->pU[x + 1];
                uint32_t i = std::upper_bound(start, end, candR[0]) - g->e1.begin();
                for (; i < g->pU[x + 1]; i++) {
                    int k = candR.idx(g->e1[i]);
                    if (k < rSize) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            }
            // printf("pU[%d]: %d ", j, pU[j]);
        }
        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu

                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    // printf("p: %d, q: %d, C[%d][%d]: %f\n", p, q, pU[u + 1] - i - 1, q - p + 1, C[pU[u + 1] - i - 1][q - p + 1]);
                    dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }
        return k;
    };
    auto buildDPTest = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);
        // todo
        // todo

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            if (rSize < g->deg1(x)) {
                for (int k = 0; k < rSize; k++) {
                    uint32_t y = candR[k];
                    if (g->connectUV(x, y)) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            } else {
                auto start = g->e1.begin() + g->pU[x];
                auto end = g->e1.begin() + g->pU[x + 1];
                uint32_t i = std::upper_bound(start, end, candR[0]) - g->e1.begin();
                for (; i < g->pU[x + 1]; i++) {
                    int k = candR.idx(g->e1[i]);
                    if (k < rSize) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            }
            // printf("pU[%d]: %d ", j, pU[j]);
        }

        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu

                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    // printf("p: %d, q: %d, C[%d][%d]: %f\n", p, q, pU[u + 1] - i - 1, q - p + 1, C[pU[u + 1] - i - 1][q - p + 1]);
                    dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }

        return std::min(k, minLR);
    };

    for (uint32_t u = 0; u < g->n1; u++) {
        if (g->deg1(u) < q) {
            continue;
        }

        uint32_t candRSize = 0;

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, candRSize);
            candRSize++;
        }

        candLtemp.clear();

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            for (uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if (w == u) {
                    // keepig the index of u so it can later get vertices larger than u
                    outPosition[i - g->pU[u]] = j;
                    break;
                }
                if (twoHopCount[w] == 0) {
                    candLtemp.push_back(w);
                }
                twoHopCount[w]++;
            }
        }

        int countCandLtemp = 0;
        for (auto w : candLtemp) {
            if (twoHopCount[w] > 1) {
                candL.changeTo(w, countCandLtemp);
                countCandLtemp++;
            }
            twoHopCount[w] = 0;
        }
        // /if no vertices found with two hop neighbors
        if (countCandLtemp == 0) {
            continue;
        }

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int candLSize = 0;
            // int outPosE = std::find(g->pV[v], g->pV[v + 1], u);
            // auto it = std::find(std::next(g->e2.begin(), g->pV[v]), std::next(g->e2.begin(), g->pV[v + 1]), u);
            // int outPosE = std::distance(g->e2.begin(), it);
            //  /assert(outPosE == outPosition[i - g->pU[u]]);
            for (uint32_t j = outPosition[i - g->pU[u]] + 1; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                // check whether the w is inside the candLTemo
                candL.changeTo(w, candLSize++);
            }

            candR.changeTo(v, --candRSize);

            if (candLSize < p - 1)
                continue;
            if (candRSize < q - 1)
                break;

            for (int i = 1; i < candRSize; i++) {
                candR.changeToByPos(i, i - 1);
            }
            // adding to shadow
            // calculate mu here
            double mu = 0;
            // printf("candLSize: %d, candRSize: %d, p: %d, q: %d\n", candLSize,
            // candRSize, p, q); sampling inside shadow to calcule mu
            int maxPossibleLength = buildDP(candLSize, candRSize, p - 1, q - 1);
            if (maxPossibleLength < minPQ) {
                continue;
            }

            double zInSubSpace = 0;
            double zInSubSpaceTest = 0;

            if (minPQ - 1 == 1) {
                for (uint32_t i = 0; i < pV[candRSize]; i++) {
                    zInSubSpace += dpV[minPQ - 1][i];
                }
            } else {
                for (uint32_t i = 0; i < pU[candLSize]; i++) {
                    zInSubSpace += dpU[minPQ - 1][i];
                }
            }
            // if (candL[0] == 555 && candR[0] == 1773) {
            //     printf("success\n");
            //     printf("zInSubSpace: %f\n", zInSubSpace);
            //     printf("candLSize: %d, candRSize: %d\n", candLSize, candRSize);
            //     for (int i = 0; i < candLSize; i++) {
            //         printf("candL[%d]: %d\n", i, candL[i]);
            //     }
            //     for (int i = 0; i < candRSize; i++) {
            //         printf("candR[%d]: %d\n", i, candR[i]);
            //     }
            //     for (int i = 0; i <= candLSize; i++) {
            //         printf("pU[%d]: %d\n", i, pU[i]);
            //     }
            //     for (int i = 0; i <= candRSize; i++) {
            //         printf("pV[%d]: %d\n", i, pV[i]);
            //     }
            //     for (int i = 0; i < pU[candLSize]; i++) {
            //         printf("e1[%d]: %d\n", i, e1[i]);
            //     }
            //     for (int i = 0; i < pV[candRSize]; i++) {
            //         printf("e2[%d]: %d\n", i, e2[i]);
            //     }
            //     exit(0);
            // }
            totalZinShadow += zInSubSpace;

            if (zInSubSpace < 1)
                continue;
            double pqBicluqeContainZ = 0;
            double auxSampleSize = 10;
            int auxSampleSize_temp = auxSampleSize;
            double t_sample = 0;

            auto auxSamplingStart = chrono::high_resolution_clock::now();
            std::pair<uint32_t, uint32_t> edge;
            std::vector<std::pair<uint32_t, uint32_t>> edges;
            vector<double> probabilities;
            vector<double> probabilitiesForFirstE;
            vector<double> Prob, probForFirstE;
            vector<uint32_t> Alias, aliasForFirstE;

            if (minPQ == 2) {
                for (int v = 0; v < candRSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ - 1][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < candLSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpU[minPQ - 1][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

            while (auxSampleSize_temp--) {
                // selecting the first edge
                std::random_device rd;

                std::mt19937 eng(rd());  // Engine (Mersenne Twister)
                std::uniform_real_distribution<> distribution(
                    0.0, 1.0);
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                // if (minPQ - 1 == 1) {
                //     for (int v = 0; v < candRSize; v++) {
                //         for (int i = pV[v]; i < pV[v + 1]; i++) {
                //             temp += dpV[minPQ - 1][mapVtoU[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(v);
                //                 selectedL.push_back(e2[i]);
                //                 preU = e2[i];
                //                 preV = v;
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // } else {
                //     double temp = 0.0;
                //     for (int u = 0; u < candLSize; u++) {
                //         for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                //             temp += dpU[minPQ - 1][mapUtoV[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(e1[i]);
                //                 selectedL.push_back(u);
                //                 preU = u;
                //                 preV = e1[i];
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // }

                int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                edge = edges[sampledIndexForFirstE];
                if (minPQ == 2) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                // printf("minPQ: %d\n", minPQ);
                // printf("---------------------\n");
                // printf("preV: %d, preU: %d\n", preV, preU);
                for (int i = 1; i < minPQ - 1; i++) {
                    // out neighborhood select u
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ - i - 1][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ - i][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }
                    // select v
                    if (i != minPQ - 2) {
                        // printf("dpV[%d][%d]: %f\n", minPQ - i - 1, mapVtoU[preE], dpV[minPQ - i - 1][mapVtoU[preE]]);
                        // printf("pU[preU]: %d, pU[preU + 1]: %d\n", pU[preU], pU[preU + 1]);
                        // printf("preV: %d\n", preV);
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ - i - 1][mapUtoV[j]];
                                // printf("temp: %f\n", temp);
                                // printf("point * dpV[minPQ - i - 1][mapVtoU[preE]]: %f\n", point * dpV[minPQ - i - 1][mapVtoU[preE]]);
                                if (temp + 1e-8 >= point * dpV[minPQ - i - 1][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (selectedL.size() != minPQ - 1) {
                    printf("selectedL.size(): %d, minPQ - 1: %d\n", selectedL.size(), minPQ - 1);
                    printf("selectedR.size(): %d, minPQ - 1: %d\n", selectedR.size(), minPQ - 1);
                    exit(0);
                }
                assert(selectedL.size() == minPQ - 1);
                int vCount = q - p + 1;
                if (minPQ - 1 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ == 2) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                // assert(pU[preU + 1] - index >= vCount);
                // assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // if (selectedR.size() != q - 1) {
                //     printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                //     printf("selectedL.size(): %d, p - 1: %d\n", selectedL.size(), p - 1);
                //     exit(0);
                // }
                assert(selectedR.size() == q - 1);
                assert(selectedL.size() == minPQ - 1);

                // check whether the sampled zstar is a biclique
                bool connected = true;
                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candL[selectedL[i]], candR[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
            }

            // calculating t_sample
            auto auxSamplingEnd = chrono::high_resolution_clock::now();
            double auxSamplingDuration = chrono::duration_cast<chrono::microseconds>(
                                             auxSamplingEnd - auxSamplingStart)
                                             .count();
            tTotal += auxSamplingDuration;
            auxTotalSamples += auxSampleSize;
            // printf("pqBicluqeContainZ: %f\n", pqBicluqeContainZ   );
            pqBicluqeEstAux += zInSubSpace * (pqBicluqeContainZ / auxSampleSize);

            // subspaceTest.push_back(Subspace(p - 1, q - 1, candL, candR, candLSize,
            // candRSize, pqBicluqeContainZ / auxSampleSize, zInSubSpace *
            // (pqBicluqeContainZ / auxSampleSize), zInSubSpace));
            std::vector<uint32_t> candLShadow;
            std::vector<uint32_t> candRShadow;
            for (int i = 0; i < candLSize; i++) {
                candLShadow.push_back(candL[i]);
            }
            for (int i = 0; i < candRSize; i++) {
                candRShadow.push_back(candR[i]);
            }

            if (zInSubSpace > 0.5) {
                shadow.push(Subspace(p - 1, q - 1, candLShadow, candRShadow, candLSize, candRSize,
                                     pqBicluqeContainZ / auxSampleSize,
                                     zInSubSpace * (pqBicluqeContainZ / auxSampleSize),
                                     zInSubSpace));
            }

            t_sample += auxSamplingDuration;

            // if (maxPossibleLength >= minPQ) {
            // 	if (pU[minPQ - 1] > 0.3) {
            // 		for (int i = 0; i < candLSize; i++) {
            // 			totalZinShadow += dpU[minPQ - 1][i];
            // 		}
            // 	}
            // }
        }
    }
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage I is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    printf("totalZinShadow: %f\n", totalZinShadow);
    auto end = chrono::high_resolution_clock::now();
    double shadowDuration =
        chrono::duration_cast<chrono::microseconds>(end - start).count();
    std::random_device rd;
    std::mt19937 eng(rd());  // Engine (Mersenne Twister)
    std::uniform_real_distribution<> distribution(
        0.0, 1.0);  // Distribution for doubles between 0.0 and 1.0 (inclusive of
    // 0.0 and exclusive of 1.0)
    vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
    printf("gamma: %f\n", gamma);
    printf("tTotal: %f\n", tTotal);
    printf("auxTotalSamples: %f\n", auxTotalSamples);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("Biclique Density: %f\n", pqBicluqeEstAux / totalZinShadow);
    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    fflush(stdout);
    printf("shadow size before the second refinement : %d \n", shadow.size());
    printf(" gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux: %f\n", gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    int count = 0;
    printf("shadowDuration: %f\n", shadowDuration);
    printf("shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux : %d\n", shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    // /shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux

    if (shadow.size() == 0) {
        printf("No biclique found fopr p: %d, q: %d, e: %f\n", p, q, e);
        return;
    }

    while (shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux) {
        Subspace min_space = shadow.top();
        if (min_space.mu > 0.5) break;
        shadow.pop();
        if (min_space.p_prime == 1 || min_space.q_prime == 1) {
            min_space.mu = 1;
            shadow.push(min_space);
        }

        uint32_t llSize = min_space.SUSize;
        uint32_t rrSize = min_space.SVSize;
        LinearSet candLL, candRR;
        candLL.resize(g->n1 + 1);
        candRR.resize(g->n2 + 1);
        for (int i = 0; i < llSize; i++) {
            candLL.changeTo(min_space.SU[i], i);
        }
        for (int i = 0; i < rrSize; i++) {
            candRR.changeTo(min_space.SV[i], i);
        }
        int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
        vector<uint32_t> pUU(g->maxDv + 1);
        vector<uint32_t> pVV(g->maxDu + 1);
        vector<uint32_t> e11(g->m);
        vector<uint32_t> e22(g->m);
        vector<uint32_t> candLtemp2(g->maxDv * g->maxDv);
        vector<uint32_t> outPosition2(g->maxDu * g->maxDu);
        std::vector<uint32_t> minSU = min_space.SU;
        std::vector<uint32_t> minSV = min_space.SV;
        std::unordered_map<uint32_t, uint32_t> minSUMap;
        std::unordered_map<uint32_t, uint32_t> minSVMap;
        vector<uint32_t> twoHopCount2(g->maxDu * g->maxDu);

        pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
        totalZinShadow -= min_space.ZCount;

        for (int i = 0; i < llSize; i++) {
            minSUMap[minSU[i]] = i;
        }
        for (int i = 0; i < rrSize; i++) {
            minSVMap[minSV[i]] = i;
        }

        /////////////////////////////
        /* ***subgraph creating*** */
        ///////////////////////////
        int edgeIndex = 0;
        for (int j = 0; j < llSize; j++) {
            uint32_t x = candLL[j];
            pUU[j + 1] = pUU[j];
            for (int k = 0; k < rrSize; k++) {
                uint32_t y = candRR[k];
                if (g->connectUV(x, y)) {
                    e11[edgeIndex] = y;
                    edgeIndex++;
                    pUU[j + 1]++;
                    pVV[k + 1]++;
                }
            }
        }

        // Resize e11 to the actual number of edges
        // Convert pVV to prefix sum array
        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }
        // Filling e22
        // edgeIndex = 0;

        int* currVV = new int[rrSize + 1];
        for (int v = 0; v <= rrSize; v++) {
            currVV[v] = pVV[v];
        }

        // Filling e22
        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                int index = minSVMap[v];
                e22[currVV[index]] = candLL[u];
                currVV[index]++;
                edgeIndex++;
            }
        }

        // No need to reset pVV since it wasn't modified
        delete[] currVV;  // Free the allocated memor

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////

        candLL.clear();
        candRR.clear();
        candLL.resize(g->n1 + 1);
        candRR.resize(g->n2 + 1);
        uint32_t totalZinThisSub = 0;
        for (uint32_t uIndex = 0; uIndex < llSize; uIndex++) {
            uint32_t u = minSU[uIndex];
            uint32_t candRSize2 = 0;

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                uint32_t v = e11[i];
                candRR.changeTo(v, candRSize2++);
            }

            candLtemp2.clear();

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                uint32_t v = e11[i];
                uint32_t vIndex = minSVMap[v];
                for (uint32_t j = pVV[vIndex + 1] - 1; j >= pVV[vIndex]; j--) {
                    uint32_t w = e22[j];
                    if (w == u) {
                        // keepig the index of u so it can later get vertices larger than u
                        outPosition2[i - pUU[uIndex]] = j;
                        break;
                    }
                    if (twoHopCount2[w] == 0) {
                        candLtemp2.push_back(w);
                    }
                    twoHopCount2[w]++;
                }
            }

            int countCandLtemp2 = 0;
            for (auto w : candLtemp2) {
                if (twoHopCount2[w] > 0) {
                    candLL.changeTo(w, countCandLtemp2++);
                }
                twoHopCount2[w] = 0;
            }

            // if no vertices found with two hop neighbors

            if (countCandLtemp2 == 0) {
                continue;
            }

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                uint32_t v = e11[i];

                int candLSize2 = 0;

                uint32_t j = outPosition2[i - pUU[uIndex]] + 1;

                uint32_t vIndex = minSVMap[v];

                for (; j < pVV[vIndex + 1]; j++) {
                    uint32_t w = e22[j];
                    // check whether the w is inside the candL
                    candLL.changeTo(w, candLSize2++);
                }

                candRR.changeTo(v, --candRSize2);

                // if (candLSize2 < min_space.p_prime - 1)
                //     continue;
                // if (candRSize2 < min_space.q_prime - 1)
                //     break;

                for (int i = 1; i < candRSize2; i++) {
                    candRR.changeToByPos(i, i - 1);
                }

                candL.clear();
                candR.clear();
                candL = candLL;
                candR = candRR;

                if (candLSize2 < min_space.p_prime - 1 || candRSize2 < min_space.q_prime - 1)
                    continue;

                int maxPossibleLength = buildDPTest(candLSize2, candRSize2, min_space.p_prime - 1, min_space.q_prime - 1);

                // if (maxPossibleLength >= minPQ2 - 1) {
                //     for (int i = 0; i < pU[candLSize2]; i++) {
                //         printf("dpU[%d][%d]: %f\n", minPQ2 - 1, i, dpU[minPQ2 - 1][i]);
                //     }
                // }
                double zInMinSubSpace = 0;
                if (minPQ2 == 2) {
                    for (uint32_t i = 0; i < pV[candRSize2]; i++) {
                        zInMinSubSpace += dpV[minPQ2 - 1][i];
                    }

                }

                else {
                    for (uint32_t i = 0; i < pU[candLSize2]; i++) {
                        zInMinSubSpace += dpU[minPQ2 - 1][i];
                    }
                }

                totalZinThisSub += zInMinSubSpace;
                totalZinShadow += zInMinSubSpace;

                if (zInMinSubSpace < 1)
                    continue;

                // sampling
                uint32_t auxSampleSize = 100;
                uint32_t auxSampleSizeTemp = auxSampleSize;
                uint32_t pqBicluqeContainZ = 0;
                std::pair<uint32_t, uint32_t> edge;
                std::vector<std::pair<uint32_t, uint32_t>> edges;
                vector<double> probabilities;
                vector<double> probabilitiesForFirstE;
                vector<double> Prob, probForFirstE;
                vector<uint32_t> Alias, aliasForFirstE;

                if (minPQ2 == 2) {
                    for (int v = 0; v < candRSize2; v++) {
                        for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                            edge.first = i;
                            edge.second = v;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpV[minPQ2 - 1][mapVtoU[i]] / zInMinSubSpace);
                        }
                    }
                } else {
                    for (int u = 0; u < candLSize2; u++) {
                        for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                            edge.first = u;
                            edge.second = i;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpU[minPQ2 - 1][mapUtoV[i]] / zInMinSubSpace);
                        }
                    }
                }

                initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

                while (auxSampleSizeTemp--) {
                    // selecting the first edge
                    selectedL.clear();
                    selectedR.clear();
                    double point = distribution(eng);
                    double temp = 0.0;
                    uint32_t preU = 0, preV = 0, preE = 0;
                    bool selectedLEmpty = true;

                    // //pasted
                    // for (int u = 0; u < candLSize2; u++) {
                    // 	for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                    // 		temp += dpU[minPQ2 - 1][mapUtoV[i]];
                    // 		if (temp + 1e-8 >= point * zInMinSubSpace) {
                    // 			selectedR.push_back(e1[i]);
                    // 			selectedL.push_back(u);
                    // 			preU = u;
                    // 			preV = e1[i];
                    // 			preE = i;
                    // 			selectedLEmpty = false;
                    // 			break;
                    // 		}
                    // 	}
                    // 	if (!selectedLEmpty) {
                    // 		break;
                    // 	}
                    // }
                    int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                    edge = edges[sampledIndexForFirstE];
                    if (minPQ2 == 2) {
                        selectedL.push_back(e2[edge.first]);
                        selectedR.push_back(edge.second);
                        preE = edge.first;
                        preU = e2[edge.first];
                        preV = edge.second;
                    } else {
                        selectedR.push_back(e1[edge.second]);
                        selectedL.push_back(edge.first);
                        preE = edge.second;
                        preU = edge.first;
                        preV = e1[edge.second];
                    }

                    assert(selectedR.size() > 0);
                    assert(selectedL.size() > 0);
                    // selecting the rest of the edge
                    for (int i = 1; i < minPQ2 - 1; i++) {
                        // out neighborhood select u
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                            if (e2[j] > preU) {
                                temp += dpV[minPQ2 - i - 1][mapVtoU[j]];
                                if (temp + 1e-8 >= point * dpU[minPQ2 - i][mapUtoV[preE]]) {
                                    uint32_t u = e2[j];
                                    selectedL.push_back(u);
                                    preU = u;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                        // select v
                        if (i != minPQ2 - 2) {
                            temp = 0.0;
                            point = distribution(eng);
                            for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                                if (e1[j] > preV) {
                                    temp += dpU[minPQ2 - i - 1][mapUtoV[j]];
                                    if (temp + 1e-8 >= point * dpV[minPQ2 - i - 1][mapVtoU[preE]]) {
                                        uint32_t v = e1[j];
                                        selectedR.push_back(v);
                                        preV = v;
                                        preE = j;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    int vCount = q - p + 1;
                    if (minPQ2 == 2) {
                        vCount = q - p;
                    }

                    uint32_t index = 0;
                    if (minPQ == 2) {
                        index = pU[preU];
                    } else {
                        auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                        index = std::distance(e1.begin(), it);
                    }

                    assert(index >= pU[preU] && index < pU[preU + 1]);

                    std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                    for (int i = 0; i < selectedRStar.size(); i++) {
                        selectedR.push_back(e1[selectedRStar[i]]);
                    }

                    assert(selectedR.size() == min_space.q_prime - 1);
                    assert(selectedL.size() == minPQ2 - 1);
                    // check whether the sampled z is a biclique
                    bool connected = true;

                    for (int i = 0; i < selectedL.size(); i++) {
                        for (int j = 0; j < selectedR.size(); j++) {
                            if (i == j)
                                continue;
                            if (i > 0 && i == j + 1)
                                continue;
                            if (!g->connectUV(candL[selectedL[i]], candR[selectedR[j]])) {
                                connected = false;
                                break;
                            }
                        }
                        if (!connected)
                            break;
                    }
                    if (!connected)
                        continue;
                    pqBicluqeContainZ++;
                }
                std::vector<uint32_t> candLShadow;
                std::vector<uint32_t> candRShadow;
                for (int i = 0; i < candLSize2; i++) {
                    candLShadow.push_back(candLL[i]);
                }
                for (int i = 0; i < candRSize2; i++) {
                    candRShadow.push_back(candRR[i]);
                }
                shadow.push(Subspace(
                    min_space.p_prime - 1, min_space.q_prime - 1, candLShadow, candRShadow,
                    candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize,
                    zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize),
                    zInMinSubSpace));

                pqBicluqeEstAux +=
                    zInMinSubSpace * pqBicluqeContainZ / (double)auxSampleSize;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        shadowDuration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        candLL.clear();
        candRR.clear();
        e11.clear();
        e22.clear();
        pVV.clear();
        outPosition2.clear();
        twoHopCount2.clear();
    }

    printf("end the second phrase\n");
    printf("shadow size: %d\n", shadow.size());

    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    printf("totalZinShadow: %f\n", totalZinShadow);
    double roughPQEstDens = pqBicluqeEstAux / totalZinShadow;
    vector<Subspace> subs;
    int batchSize = gamma / roughPQEstDens;
    std::pair<uint32_t, uint32_t> edge;
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    printf("roughPQEstDens: %f\n", roughPQEstDens);
    vector<double> probabilities;
    vector<double> probabilitiesForFirstE;
    vector<Subspace> shadowVec;
    vector<double> Prob, probForFirstE;
    vector<uint32_t> Alias, aliasForFirstE;
    std::mt19937 gen(rd());
    uint32_t pqBicluqeContainZ = 0;

    printf("shadow.size(): %d\n", shadow.size());

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage II is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////

    while (!shadow.empty()) {
        shadowVec.push_back(shadow.top());
        shadow.pop();
    }

    printf("totalZinShadow: %f\n", totalZinShadow);

    for (const auto& subspace : shadowVec) {
        probabilities.push_back(static_cast<double>(subspace.ZCount) / totalZinShadow);
    }

    initializeAliasMethod(probabilities, Prob, Alias);
    int sampledIndex = getIndexAlias(gen, Prob, Alias);
    std::unordered_map<int, int> subSpaceAndC;
    printf("batchSize: %d\n", batchSize);
    int successSamples = 0;
    int t = 0;
    bool gammaReached = false;
    double totalZFinal = 0;

    while (successSamples < int(gamma)) {
        subSpaceAndC.clear();

        for (int i = 0; i < batchSize; i++) {
            int sampledIndex = getIndexAlias(gen, Prob, Alias);
            if (subSpaceAndC.count(sampledIndex) == 0) {
                subSpaceAndC[sampledIndex] = 1;
            } else {
                subSpaceAndC[sampledIndex] = subSpaceAndC[sampledIndex] + 1;
            }
        }

        for (const auto& batch : subSpaceAndC) {
            Subspace s = shadowVec[batch.first];
            totalZFinal += s.ZCount;
            candL.clear();
            candR.clear();
            candL.resize(g->n1 + 1);
            candR.resize(g->n2 + 1);
            // candL = s.SU;
            // candR = s.SV;
            // printf("s.SU.size(): %d\n", s.SU.size());
            // printf("s.SUSize: %d\n", s.SUSize);
            // for (int i = 0; i < s.SUSize; i++) {
            //     printf("s.SU[i]: %d\n", s.SU[i]);
            // }
            // continue;
            // exit(0);
            for (int i = 0; i < s.SUSize; i++) {
                candL.changeTo(s.SU[i], i);
            }
            for (int i = 0; i < s.SVSize; i++) {
                candR.changeTo(s.SV[i], i);
            }

            int maxLength = buildDP(s.SUSize, s.SVSize, s.p_prime, s.q_prime);
            uint32_t llSize = s.SUSize;
            uint32_t rrSize = s.SVSize;
            LinearSet candLL, candRR;
            candLL.resize(g->n1 + 1);
            candRR.resize(g->n2 + 1);
            // candLL = s.SU;
            // candRR = s.SV;
            for (int i = 0; i < s.SUSize; i++) {
                candLL.changeTo(s.SU[i], i);
            }
            for (int i = 0; i < s.SVSize; i++) {
                candRR.changeTo(s.SV[i], i);
            }
            int minPQ2 = std::min(s.p_prime, s.q_prime);
            double zInSubSpace = s.ZCount;

            edges.clear();
            probabilitiesForFirstE.clear();
            probForFirstE.clear();
            aliasForFirstE.clear();

            if (minPQ2 == 1) {
                for (int v = 0; v < s.SVSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ2][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < s.SUSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        // printf("dpU[minPQ2][mapUtoV[i]]: %f\n", dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                        probabilitiesForFirstE.push_back(dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }
            double sum = 0;
            for (double d : probabilitiesForFirstE) {
                sum += d;
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);
            t += batch.second;
            for (int i = 0; i < batch.second; i++) {
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                int sampledIndexForFirstE = getIndexAlias(gen, probForFirstE, aliasForFirstE);
                edge = edges[sampledIndexForFirstE];

                if (minPQ2 == 1) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                for (int i = 1; i < minPQ2; i++) {
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ2 - i][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ2 - i + 1][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }

                    // select v
                    if (i != minPQ2 - 1) {
                        temp = 0.0;
                        point = distribution(eng);
                        bool f = false;
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ2 - i][mapUtoV[j]];
                                if (temp + 1e-8 >= point * dpV[minPQ2 - i][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    f = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                int vCount = q - p + 1;
                if (minPQ2 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ2 == 1) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // for (size_t i = 0; i < selectedL.size(); i++) {
                //     printf("%d ", selectedL[i]);
                // }
                // printf("\n");
                assert(selectedR.size() == s.q_prime);
                assert(selectedL.size() == minPQ2);
                // check whether the sample is a biclique
                bool connected = true;

                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candL[selectedL[i]], candR[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
                successSamples++;

                if (successSamples >= (int)gamma) {
                    gammaReached = true;
                    break;
                }
            }

            if (gammaReached) break;
        }
    }

    printf("pqBicluqeContainZ: %d\n", pqBicluqeContainZ);
    printf("t: %d\n", t);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("successSamples: %d\n", successSamples);
    printf("density of biclique: %f\n", pqBicluqeContainZ / (double)t);
    printf("biclique count: %f\n", pqBicluqeContainZ / (double)t * totalZinShadow);
    printf("time: %f\n", (double)(std::chrono::high_resolution_clock::now() - start).count() / 1e9);
}

// Shadow Builder with vectors
void accuracy::shadowBuilderZStar2(int p, int q, double epsilon) {
    std::signal(SIGSEGV, segfault_handler);
    if (p > q) {
        swap(p, q);
    }
    printf(" Working on p q e: %d - %d - %.2f \n", p, q, epsilon);
    auto start = chrono::high_resolution_clock::now();
    minPQ = min(p, q);

    // get the total structure count in the input graph
    // e1 hols u's neighbors

    std::vector<uint32_t> candLtemp(g->maxDv);
    std::vector<uint32_t> outPosition(g->maxDu);
    vector<int> twoHopCount(g->n1);
    candL.resize(g->n1);
    candR.resize(g->n2);
    std::vector<uint32_t> candLVec, candRVec;
    double totalZinShadow = 0;
    double totalZinShadowTest = 0;

    uint32_t maxE = std::min((1u << 21), g->maxDu * g->maxDv);
    // for Shadow
    uint32_t cnt_PQ = 1;

    priority_queue<Subspace, std::vector<Subspace>, CompareSubspace> shadow;
    vector<Subspace> subspaceTest;
    uint32_t totalS = 0;

    std::vector<uint32_t> pV(g->maxDu);  // Adjusted size to accommodate pV[v + 1]
    std::vector<uint32_t> pU(g->maxDv);
    std::vector<uint32_t> e1(maxE);  // Assuming maximum possible size
    std::vector<uint32_t> e2(
        maxE);  // Assuming maximum possible sizmaxPossibleLengthe
    std::vector<double> ddp(maxE, 0.0);
    std::vector<std::vector<double>> dpU(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<std::vector<double>> dpV(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<uint32_t> mapUtoV(maxE);
    std::vector<uint32_t> mapVtoU(maxE);
    const double delta = 0.05;
    double gamma = 4 * (1 + epsilon) * (std::exp(1) - 2) * std::log(2 / delta) /
                   (epsilon * epsilon);

    double tSample = 0;
    double tTotal = 0;
    double auxTotalSamples = 1000;
    double pqBicluqeEstAux = 1;
    int candLSum = 0;
    int candRSum = 0;

    for (uint32_t e = 0; e < maxE; e++) {
        dpU[0][e] = 0;
        dpU[1][e] = 1;
        dpV[0][e] = 1;
    }

    // DP Buidling
    auto buildDP = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            if (rSize < g->deg1(x)) {
                for (int k = 0; k < rSize; k++) {
                    uint32_t y = candR[k];
                    if (g->connectUV(x, y)) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            } else {
                auto start = g->e1.begin() + g->pU[x];
                auto end = g->e1.begin() + g->pU[x + 1];
                uint32_t i = std::upper_bound(start, end, candR[0]) - g->e1.begin();
                for (; i < g->pU[x + 1]; i++) {
                    int k = candR.idx(g->e1[i]);
                    if (k < rSize) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            }
            // printf("pU[%d]: %d ", j, pU[j]);
        }
        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu

                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    // printf("p: %d, q: %d, C[%d][%d]: %f\n", p, q, pU[u + 1] - i - 1, q - p + 1, C[pU[u + 1] - i - 1][q - p + 1]);
                    dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }
        return k;
    };
    auto buildDP2Vec = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candLVec[j];
            pU[j + 1] = pU[j];
            for (int k = 0; k < rSize; k++) {
                uint32_t y = candRVec[k];
                if (g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
        }

        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }
        for (uint32_t u = 0; u < lSize; u++) {
            for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu
                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    seg_value = dpV[1].size();
                    pU_seg = i;
                    pU1_seg = maxE;
                    if (pU[u + 1] - i - 1 == 0) {
                        dpV[1][i] = 0;
                    } else {
                        double y = C[pU[u + 1] - i - 1][q - p + 1];
                        dpV[1][i] = y;
                    }
                    // dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }

        return std::min(k, minLR);
    };

    for (uint32_t u = 0; u < g->n1; u++) {
        if (g->deg1(u) < q) {
            continue;
        }

        uint32_t candRSize = 0;

        // for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
        //     uint32_t v = g->e1[i];
        //     candR.changeTo(v, candRSize);
        //     candRSize++;
        // }

        candLtemp.clear();

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, candRSize);
            candRSize++;
            for (uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if (w == u) {
                    // keepig the index of u so it can later get vertices larger than u
                    outPosition[i - g->pU[u]] = j;
                    break;
                }
                if (twoHopCount[w] == 0) {
                    candLtemp.push_back(w);
                }
                twoHopCount[w]++;
            }
        }

        int countCandLtemp = 0;
        for (auto w : candLtemp) {
            if (twoHopCount[w] > 1) {
                candL.changeTo(w, countCandLtemp);
                countCandLtemp++;
            }
            twoHopCount[w] = 0;
        }
        // /if no vertices found with two hop neighbors
        if (countCandLtemp == 0) {
            continue;
        }

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            assert(v == candR[0]);
            // if (v != candR[0]) {
            //     printf("v: %d\n", v);
            //     for (int i = 0; i < candRSize; i++) {
            //         printf("candR[%d]: %d\n", i, candR[i]);
            //     }
            // }
            int candLSize = 0;
            // int outPosE = std::find(g->pV[v], g->pV[v + 1], u);
            // auto it = std::find(std::next(g->e2.begin(), g->pV[v]), std::next(g->e2.begin(), g->pV[v + 1]), u);
            // int outPosE = std::distance(g->e2.begin(), it);
            //  /assert(outPosE == outPosition[i - g->pU[u]]);
            for (uint32_t j = outPosition[i - g->pU[u]] + 1; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                // check whether the w is inside the candLTemp
                if (candL.idx(w) < countCandLtemp) candL.changeTo(w, candLSize++);
            }

            candR.changeTo(v, --candRSize);

            // if (candLSize < p - 1)
            //     continue;
            // if (candRSize < q - 1)
            //     break;
            for (int i = 1; i < candRSize; i++) {
                candR.changeToByPos(i, i - 1);
            }
            if (candRSize < q - 1)
                break;
            if (candLSize < p - 1)
                continue;

            // adding to shadow
            // calculate mu here
            double mu = 0;
            // printf("candLSize: %d, candRSize: %d, p: %d, q: %d\n", candLSize,
            // candRSize, p, q); sampling inside shadow to calcule mu
            int maxPossibleLength = buildDP(candLSize, candRSize, p - 1, q - 1);
            if (maxPossibleLength < minPQ) {
                continue;
            }

            double zInSubSpace = 0;
            double zInSubSpaceTest = 0;

            if (minPQ - 1 == 1) {
                for (uint32_t i = 0; i < pV[candRSize]; i++) {
                    zInSubSpace += dpV[minPQ - 1][i];
                }
            } else {
                for (uint32_t i = 0; i < pU[candLSize]; i++) {
                    zInSubSpace += dpU[minPQ - 1][i];
                }
            }
            // if (candL[0] == 555 && candR[0] == 1773) {
            //     printf("success\n");
            //     printf("zInSubSpace: %f\n", zInSubSpace);
            //     printf("candLSize: %d, candRSize: %d\n", candLSize, candRSize);
            //     for (int i = 0; i < candLSize; i++) {
            //         printf("candL[%d]: %d\n", i, candL[i]);
            //     }
            //     for (int i = 0; i < candRSize; i++) {
            //         printf("candR[%d]: %d\n", i, candR[i]);
            //     }
            //     for (int i = 0; i <= candLSize; i++) {
            //         printf("pU[%d]: %d\n", i, pU[i]);
            //     }
            //     for (int i = 0; i <= candRSize; i++) {
            //         printf("pV[%d]: %d\n", i, pV[i]);
            //     }
            //     for (int i = 0; i < pU[candLSize]; i++) {
            //         printf("e1[%d]: %d\n", i, e1[i]);
            //     }
            //     for (int i = 0; i < pV[candRSize]; i++) {
            //         printf("e2[%d]: %d\n", i, e2[i]);
            //     }
            //     exit(0);
            // }
            totalZinShadow += zInSubSpace;

            if (zInSubSpace < 1)
                continue;
            double pqBicluqeContainZ = 0;
            double auxSampleSize = 10;
            int auxSampleSize_temp = auxSampleSize;
            double t_sample = 0;

            auto auxSamplingStart = chrono::high_resolution_clock::now();
            std::pair<uint32_t, uint32_t> edge;
            std::vector<std::pair<uint32_t, uint32_t>> edges;
            vector<double> probabilities;
            vector<double> probabilitiesForFirstE;
            vector<double> Prob, probForFirstE;
            vector<uint32_t> Alias, aliasForFirstE;

            if (minPQ == 2) {
                for (int v = 0; v < candRSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ - 1][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < candLSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpU[minPQ - 1][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

            while (auxSampleSize_temp--) {
                // selecting the first edge
                std::random_device rd;

                std::mt19937 eng(rd());  // Engine (Mersenne Twister)
                std::uniform_real_distribution<> distribution(
                    0.0, 1.0);
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                // if (minPQ - 1 == 1) {
                //     for (int v = 0; v < candRSize; v++) {
                //         for (int i = pV[v]; i < pV[v + 1]; i++) {
                //             temp += dpV[minPQ - 1][mapVtoU[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(v);
                //                 selectedL.push_back(e2[i]);
                //                 preU = e2[i];
                //                 preV = v;
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // } else {
                //     double temp = 0.0;
                //     for (int u = 0; u < candLSize; u++) {
                //         for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                //             temp += dpU[minPQ - 1][mapUtoV[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(e1[i]);
                //                 selectedL.push_back(u);
                //                 preU = u;
                //                 preV = e1[i];
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // }

                int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                edge = edges[sampledIndexForFirstE];
                if (minPQ == 2) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                // printf("minPQ: %d\n", minPQ);
                // printf("---------------------\n");
                // printf("preV: %d, preU: %d\n", preV, preU);
                for (int i = 1; i < minPQ - 1; i++) {
                    // out neighborhood select u
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ - i - 1][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ - i][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }
                    // select v
                    if (i != minPQ - 2) {
                        // printf("dpV[%d][%d]: %f\n", minPQ - i - 1, mapVtoU[preE], dpV[minPQ - i - 1][mapVtoU[preE]]);
                        // printf("pU[preU]: %d, pU[preU + 1]: %d\n", pU[preU], pU[preU + 1]);
                        // printf("preV: %d\n", preV);
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ - i - 1][mapUtoV[j]];
                                // printf("temp: %f\n", temp);
                                // printf("point * dpV[minPQ - i - 1][mapVtoU[preE]]: %f\n", point * dpV[minPQ - i - 1][mapVtoU[preE]]);
                                if (temp + 1e-8 >= point * dpV[minPQ - i - 1][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (selectedL.size() != minPQ - 1) {
                    printf("selectedL.size(): %d, minPQ - 1: %d\n", selectedL.size(), minPQ - 1);
                    continue;
                }
                assert(selectedL.size() == minPQ - 1);
                int vCount = q - p + 1;
                if (minPQ - 1 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ == 2) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                // assert(pU[preU + 1] - index >= vCount);
                // assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // if (selectedR.size() != q - 1) {
                //     printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                //     printf("selectedL.size(): %d, p - 1: %d\n", selectedL.size(), p - 1);
                //     exit(0);
                // }
                if (selectedR.size() != q - 1) {
                    printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                    continue;
                }
                assert(selectedR.size() == q - 1);
                assert(selectedL.size() == minPQ - 1);

                // check whether the sampled zstar is a biclique
                bool connected = true;
                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candL[selectedL[i]], candR[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
            }

            // calculating t_sample
            auto auxSamplingEnd = chrono::high_resolution_clock::now();
            double auxSamplingDuration = chrono::duration_cast<chrono::microseconds>(
                                             auxSamplingEnd - auxSamplingStart)
                                             .count();
            tTotal += auxSamplingDuration;
            auxTotalSamples += auxSampleSize;
            // printf("pqBicluqeContainZ: %f\n", pqBicluqeContainZ   );
            pqBicluqeEstAux += zInSubSpace * (pqBicluqeContainZ / auxSampleSize);

            // subspaceTest.push_back(Subspace(p - 1, q - 1, candL, candR, candLSize,
            // candRSize, pqBicluqeContainZ / auxSampleSize, zInSubSpace *
            // (pqBicluqeContainZ / auxSampleSize), zInSubSpace));
            std::vector<uint32_t> candLShadow;
            std::vector<uint32_t> candRShadow;
            for (int i = 0; i < candLSize; i++) {
                candLShadow.push_back(candL[i]);
            }
            for (int i = 0; i < candRSize; i++) {
                candRShadow.push_back(candR[i]);
            }

            if (zInSubSpace > 0.5) {
                shadow.push(Subspace(p - 1, q - 1, candLShadow, candRShadow, candLSize, candRSize,
                                     pqBicluqeContainZ / auxSampleSize,
                                     zInSubSpace * (pqBicluqeContainZ / auxSampleSize),
                                     zInSubSpace));
            }

            t_sample += auxSamplingDuration;

            // if (maxPossibleLength >= minPQ) {
            // 	if (pU[minPQ - 1] > 0.3) {
            // 		for (int i = 0; i < candLSize; i++) {
            // 			totalZinShadow += dpU[minPQ - 1][i];
            // 		}
            // 	}
            // }
        }
    }

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage I is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    printf("totalZinShadow: %f\n", totalZinShadow);
    auto end = chrono::high_resolution_clock::now();
    double shadowDuration =
        chrono::duration_cast<chrono::microseconds>(end - start).count();
    std::random_device rd;
    std::mt19937 eng(rd());  // Engine (Mersenne Twister)
    std::uniform_real_distribution<> distribution(
        0.0, 1.0);  // Distribution for doubles between 0.0 and 1.0 (inclusive of
    // 0.0 and exclusive of 1.0)
    vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
    printf("gamma: %f\n", gamma);
    printf("tTotal: %f\n", tTotal);
    printf("auxTotalSamples: %f\n", auxTotalSamples);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("Biclique Density: %f\n", pqBicluqeEstAux / totalZinShadow);
    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    fflush(stdout);
    printf("shadow size before the second refinement : %d \n", shadow.size());
    printf(" gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux: %f\n", gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    int count = 0;
    printf("shadowDuration: %f\n", shadowDuration);
    printf("shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux : %d\n", shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    // /shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux
    fflush(stdout);
    if (shadow.size() == 0) {
        printf("No biclique found fopr p: %d, q: %d\n", p, q);
        return;
    }

    while (shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux) {
        Subspace min_space = shadow.top();
        if (min_space.mu > 0.5) break;
        shadow.pop();
        if (min_space.p_prime == 1 || min_space.q_prime == 1) {
            min_space.mu = 1;
            shadow.push(min_space);
        }

        uint32_t llSize = min_space.SUSize;
        uint32_t rrSize = min_space.SVSize;
        std::vector<uint32_t> candLL, candRR;
        candLL = min_space.SU;
        candRR = min_space.SV;
        int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
        vector<uint32_t> pUU(g->maxDv + 1);
        vector<uint32_t> pVV(g->maxDu + 1);
        vector<uint32_t> e11(g->m);
        vector<uint32_t> e22(g->m);
        // vector<uint32_t> candLtemp2(g->maxDv * g->maxDv);
        // vector<uint32_t> outPosition2(g->maxDu * g->maxDu);
        std::vector<uint32_t> minSU = min_space.SU;
        std::vector<uint32_t> minSV = min_space.SV;
        std::unordered_map<uint32_t, uint32_t> minSUMap;
        std::unordered_map<uint32_t, uint32_t> minSVMap;
        // vector<uint32_t> twoHopCount2(g->maxDu * g->maxDu);

        pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
        totalZinShadow -= min_space.ZCount;

        for (int i = 0; i < llSize; i++) {
            minSUMap[minSU[i]] = i;
        }
        for (int i = 0; i < rrSize; i++) {
            minSVMap[minSV[i]] = i;
        }

        /////////////////////////////
        /* ***subgraph creating*** */
        ///////////////////////////
        int edgeIndex = 0;
        for (int j = 0; j < llSize; j++) {
            uint32_t x = candLL[j];
            pUU[j + 1] = pUU[j];
            for (int k = 0; k < rrSize; k++) {
                uint32_t y = candRR[k];
                if (g->connectUV(x, y)) {
                    e11[edgeIndex++] = y;
                    pUU[j + 1]++;
                    pVV[k + 1]++;
                }
            }
        }

        // Resize e11 to the actual number of edges
        // Convert pVV to prefix sum array
        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }
        assert(pVV[rrSize] == pUU[llSize]);
        // for (int i = 0; i <= llSize; i++) {
        //     printf("pUU[%d]: %d\n", i, pUU[i]);
        // }
        // for (int i = 0; i <= rrSize; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        // Filling e22
        // edgeIndex = 0;

        // Filling e22
        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                int index = minSVMap[v];
                e22[pVV[index]] = candLL[u];
                pVV[index]++;
            }
        }

        // Reset pVV to prefix sum array
        for (int i = rrSize + 1; i >= 1; i--) {
            pVV[i] = pVV[i - 1];
        }
        pVV[0] = 0;
        // if (pUU[llSize] != pVV[rrSize]) {
        //     printf("minSpace.p_prime: %d, minSpace.q_prime: %d\n", min_space.p_prime, min_space.q_prime);
        //     printf("--------------------\n");

        //     for (int i = 0; i <= llSize; i++) {
        //         printf("pUU[%d]: %d\n", i, pUU[i]);
        //     }
        //     for (int i = 0; i <= rrSize; i++) {
        //         printf("pVV[%d]: %d\n", i, pVV[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < pUU[llSize]; i++) {
        //         printf("e11[%d]: %d\n", i, e11[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < pVV[rrSize]; i++) {
        //         printf("e22[%d]: %d\n", i, e22[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < llSize; i++) {
        //         printf("candLL[%d]: %d\n", i, candLL[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < rrSize; i++) {
        //         printf("candRR[%d]: %d\n", i, candRR[i]);
        //     }
        //     exit(0);
        // }
        assert(pVV[rrSize] == pUU[llSize]);

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////
        // printf("min_space.subgraph: %d\n", min_space.subgraph);
        // printf("min_space.p_prime: %d\n", min_space.p_prime);
        // for (int i = 0; i < llSize; i++) {
        //     printf("candLL[%d]: %d\n", i, candLL[i]);
        // }
        // for (int i = 0; i < rrSize; i++) {
        //     printf("candRR[%d]: %d\n", i, candRR[i]);
        // }
        // for (int i = 0; i <= llSize; i++) {
        //     printf("pUU[%d]: %d\n", i, pUU[i]);
        // }
        // for (int i = 0; i <= rrSize; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        // for (int i = 0; i <= pUU[llSize]; i++) {
        //     printf("e11[%d]: %d\n", i, e11[i]);
        // }
        // for (int i = 0; i <= pVV[rrSize]; i++) {
        //     printf("e22[%d]: %d\n", i, e22[i]);
        // }
        candLL.clear();
        candRR.clear();

        uint32_t totalZinThisSub = 0;
        for (uint32_t uIndex = 0; uIndex < llSize; uIndex++) {
            candRR.clear();
            uint32_t u = minSU[uIndex];
            uint32_t candRSize2 = 0;

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                uint32_t v = e11[i];
                candRR.push_back(v);
                candRSize2++;
            }
            // for (int i = 0; i < candRSize2; i++) {
            //     printf(" main loop candRR[%d]: %d\n", i, candRR[i]);
            // }
            // candLtemp2.clear();

            // for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
            //     uint32_t v = e11[i];
            //     uint32_t vIndex = minSVMap[v];
            //     for (uint32_t j = pVV[vIndex + 1] - 1; j >= pVV[vIndex]; j--) {
            //         uint32_t w = e22[j];
            //         if (w == u) {
            //             // keepig the index of u so it can later get vertices larger than u
            //             outPosition2[i - pUU[uIndex]] = j;
            //             break;
            //         }
            //         // if (twoHopCount2[w] == 0) {
            //         //     candLtemp2.push_back(w);
            //         // }
            //         // twoHopCount2[w]++;
            //     }
            // }

            // int countCandLtemp2 = 0;
            // for (auto w : candLtemp2) {
            //     if (twoHopCount2[w] > 0) {
            //         candLL.push_back(w);
            //         countCandLtemp2++;
            //     }
            //     twoHopCount2[w] = 0;
            // }

            // // if no vertices found with two hop neighbors

            // if (countCandLtemp2 == 0) {
            //     continue;
            // }

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                candLL.clear();
                uint32_t v = e11[i];
                assert(v == candRR[0]);
                int candLSize2 = 0;

                // uint32_t j = outPosition2[i - pUU[uIndex]] + 1;
                uint32_t j = std::find(e22.begin() + pVV[minSVMap[v]], e22.begin() + pVV[minSVMap[v] + 1], u) - e22.begin() + 1;
                uint32_t vIndex = minSVMap[v];

                for (; j < pVV[vIndex + 1]; j++) {
                    uint32_t w = e22[j];
                    // check whether the w is inside the candL
                    candLL.push_back(w);
                    candLSize2++;
                }
                // fflush(stdout);
                // if (candRR[0] != v) {
                //     printf("candR is not equal to v\n");
                //     printf("v: %d\n", v);
                //     for (int i = 0; i < candRSize2; i++) {
                //         printf("candRR[%d]: %d\n", i, candRR[i]);
                //     }
                //     fflush(stdout);
                //     exit(0);
                // }

                //  candRR.changeTo(v, --candRSize2);
                // printf("candR befdore: %d\n", candRSize2);
                // for (int i = 0; i < candRSize2; i++) {
                //     printf("candRR[%d]: %d\n", i, candRR[i]);
                // }
                candRSize2--;
                candRR.erase(candRR.begin());
                // printf("candR after: %d\n", candRSize2);
                // for (int i = 0; i < candRSize2; i++) {
                //     printf("candRR[%d]: %d\n", i, candRR[i]);
                // }
                if (candRSize2 < min_space.q_prime - 1) {
                    break;
                }
                if (candLSize2 < min_space.p_prime - 1) {
                    continue;
                }

                // for (int i = 0; i < candRSize2; i++) {
                //     candRR.changeToByPos(i, i - 1);
                // }

                candLVec = candLL;
                candRVec = candRR;

                int maxPossibleLength = buildDP2Vec(candLSize2, candRSize2, min_space.p_prime - 1, min_space.q_prime - 1);

                double zInMinSubSpace = 0;
                if (minPQ2 == 2) {
                    for (uint32_t i = 0; i < pV[candRSize2]; i++) {
                        zInMinSubSpace += dpV[minPQ2 - 1][i];
                    }

                }

                else {
                    for (uint32_t i = 0; i < pU[candLSize2]; i++) {
                        zInMinSubSpace += dpU[minPQ2 - 1][i];
                    }
                }

                totalZinThisSub += zInMinSubSpace;
                totalZinShadow += zInMinSubSpace;

                if (zInMinSubSpace < 1)
                    continue;

                // sampling
                uint32_t auxSampleSize = 100;
                uint32_t auxSampleSizeTemp = auxSampleSize;
                uint32_t pqBicluqeContainZ = 0;
                std::pair<uint32_t, uint32_t> edge;
                std::vector<std::pair<uint32_t, uint32_t>> edges;
                vector<double> probabilities;
                vector<double> probabilitiesForFirstE;
                vector<double> Prob, probForFirstE;
                vector<uint32_t> Alias, aliasForFirstE;

                if (minPQ2 == 2) {
                    for (int v = 0; v < candRSize2; v++) {
                        for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                            edge.first = i;
                            edge.second = v;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpV[minPQ2 - 1][mapVtoU[i]] / zInMinSubSpace);
                        }
                    }
                } else {
                    for (int u = 0; u < candLSize2; u++) {
                        for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                            edge.first = u;
                            edge.second = i;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpU[minPQ2 - 1][mapUtoV[i]] / zInMinSubSpace);
                        }
                    }
                }

                initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

                while (auxSampleSizeTemp--) {
                    // selecting the first edge
                    selectedL.clear();
                    selectedR.clear();
                    double point = distribution(eng);
                    double temp = 0.0;
                    uint32_t preU = 0, preV = 0, preE = 0;
                    bool selectedLEmpty = true;

                    // //pasted
                    // for (int u = 0; u < candLSize2; u++) {
                    // 	for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                    // 		temp += dpU[minPQ2 - 1][mapUtoV[i]];
                    // 		if (temp + 1e-8 >= point * zInMinSubSpace) {
                    // 			selectedR.push_back(e1[i]);
                    // 			selectedL.push_back(u);
                    // 			preU = u;
                    // 			preV = e1[i];
                    // 			preE = i;
                    // 			selectedLEmpty = false;
                    // 			break;
                    // 		}
                    // 	}
                    // 	if (!selectedLEmpty) {
                    // 		break;
                    // 	}
                    // }
                    int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                    edge = edges[sampledIndexForFirstE];
                    if (minPQ2 == 2) {
                        selectedL.push_back(e2[edge.first]);
                        selectedR.push_back(edge.second);
                        preE = edge.first;
                        preU = e2[edge.first];
                        preV = edge.second;
                    } else {
                        selectedR.push_back(e1[edge.second]);
                        selectedL.push_back(edge.first);
                        preE = edge.second;
                        preU = edge.first;
                        preV = e1[edge.second];
                    }

                    assert(selectedR.size() > 0);
                    assert(selectedL.size() > 0);
                    // selecting the rest of the edge
                    for (int i = 1; i < minPQ2 - 1; i++) {
                        // out neighborhood select u
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                            if (e2[j] > preU) {
                                temp += dpV[minPQ2 - i - 1][mapVtoU[j]];
                                if (temp + 1e-8 >= point * dpU[minPQ2 - i][mapUtoV[preE]]) {
                                    uint32_t u = e2[j];
                                    selectedL.push_back(u);
                                    preU = u;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                        // select v
                        if (i != minPQ2 - 2) {
                            temp = 0.0;
                            point = distribution(eng);
                            for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                                if (e1[j] > preV) {
                                    temp += dpU[minPQ2 - i - 1][mapUtoV[j]];
                                    if (temp + 1e-8 >= point * dpV[minPQ2 - i - 1][mapVtoU[preE]]) {
                                        uint32_t v = e1[j];
                                        selectedR.push_back(v);
                                        preV = v;
                                        preE = j;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    int vCount = q - p + 1;
                    if (minPQ2 == 2) {
                        vCount = q - p;
                    }

                    uint32_t index = 0;
                    if (minPQ == 2) {
                        index = pU[preU];
                    } else {
                        auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                        index = std::distance(e1.begin(), it);
                    }

                    assert(index >= pU[preU] && index < pU[preU + 1]);

                    std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                    for (int i = 0; i < selectedRStar.size(); i++) {
                        selectedR.push_back(e1[selectedRStar[i]]);
                    }
                    if (selectedR.size() != min_space.q_prime - 1 || selectedL.size() != minPQ2 - 1) {
                        printf("issue is sample zsatar sizes \n");
                        printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                        continue;
                    }
                    assert(selectedR.size() == min_space.q_prime - 1);
                    assert(selectedL.size() == minPQ2 - 1);
                    // check whether the sampled z is a biclique
                    bool connected = true;

                    for (int i = 0; i < selectedL.size(); i++) {
                        for (int j = 0; j < selectedR.size(); j++) {
                            if (i == j)
                                continue;
                            if (i > 0 && i == j + 1)
                                continue;
                            if (!g->connectUV(candLVec[selectedL[i]], candRVec[selectedR[j]])) {
                                connected = false;
                                break;
                            }
                        }
                        if (!connected)
                            break;
                    }
                    if (!connected)
                        continue;
                    pqBicluqeContainZ++;
                }
                Subspace s = Subspace(min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR, candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize, zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize), zInMinSubSpace);
                s.subgraph = s.subgraph + 1;
                shadow.push(s);
                // shadow.push(Subspace(
                //     min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR,
                //     candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize,
                //     zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize),
                //     zInMinSubSpace));

                pqBicluqeEstAux +=
                    zInMinSubSpace * pqBicluqeContainZ / (double)auxSampleSize;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        shadowDuration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        candLL.clear();
        candRR.clear();
        e11.clear();
        e22.clear();
        pVV.clear();
        // outPosition2.clear();
        // twoHopCount2.clear();
    }

    printf("end the second phrase\n");
    printf("shadow size: %d\n", shadow.size());

    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    printf("totalZinShadow: %f\n", totalZinShadow);
    double roughPQEstDens = pqBicluqeEstAux / totalZinShadow;
    vector<Subspace> subs;
    uint32_t batchSize = gamma / roughPQEstDens;
    std::pair<uint32_t, uint32_t> edge;
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    printf("roughPQEstDens: %f\n", roughPQEstDens);
    vector<double> probabilities;
    vector<double> probabilitiesForFirstE;
    vector<Subspace> shadowVec;
    vector<double> Prob, probForFirstE;
    vector<uint32_t> Alias, aliasForFirstE;
    std::mt19937 gen(rd());
    uint32_t pqBicluqeContainZ = 0;

    printf("shadow.size(): %d\n", shadow.size());
    if (shadow.size() == 0) {
        printf("No biclique found fopr p: %d, q: %d\n", p, q);
        return;
    }
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage II is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////

    while (!shadow.empty()) {
        shadowVec.push_back(shadow.top());
        shadow.pop();
    }

    printf("totalZinShadow: %f\n", totalZinShadow);

    for (const auto& subspace : shadowVec) {
        probabilities.push_back(static_cast<double>(subspace.ZCount) / totalZinShadow);
    }

    initializeAliasMethod(probabilities, Prob, Alias);
    int sampledIndex = getIndexAlias(gen, Prob, Alias);
    std::unordered_map<int, int> subSpaceAndC;
    printf("batchSize: %d\n", batchSize);
    int successSamples = 0;
    int t = 0;
    bool gammaReached = false;
    double totalZFinal = 0;

    while (successSamples < int(gamma)) {
        subSpaceAndC.clear();

        for (int i = 0; i < batchSize; i++) {
            int sampledIndex = getIndexAlias(gen, Prob, Alias);
            if (subSpaceAndC.count(sampledIndex) == 0) {
                subSpaceAndC[sampledIndex] = 1;
            } else {
                subSpaceAndC[sampledIndex] = subSpaceAndC[sampledIndex] + 1;
            }
        }

        for (const auto& batch : subSpaceAndC) {
            Subspace s = shadowVec[batch.first];
            totalZFinal += s.ZCount;
            candLVec.clear();
            candRVec.clear();
            // candL.resize(g->n1 + 1);
            // candR.resize(g->n2 + 1);
            candLVec = s.SU;
            candRVec = s.SV;
            // printf("s.SU.size(): %d\n", s.SU.size());
            // printf("s.SUSize: %d\n", s.SUSize);
            // for (int i = 0; i < s.SUSize; i++) {
            //     printf("s.SU[i]: %d\n", s.SU[i]);
            // }
            // continue;
            // exit(0);
            // for (int i = 0; i < s.SUSize; i++) {
            //     candL.changeTo(s.SU[i], i);
            // }
            // for (int i = 0; i < s.SVSize; i++) {
            //     candR.changeTo(s.SV[i], i);
            // }

            int maxLength = buildDP2Vec(s.SUSize, s.SVSize, s.p_prime, s.q_prime);
            uint32_t llSize = s.SUSize;
            uint32_t rrSize = s.SVSize;
            // LinearSet candLL, candRR;
            // candLL.resize(g->n1 + 1);
            // candRR.resize(g->n2 + 1);
            // candLL = s.SU;
            // candRR = s.SV;
            // for (int i = 0; i < s.SUSize; i++) {
            //     candLL.changeTo(s.SU[i], i);
            // }
            // for (int i = 0; i < s.SVSize; i++) {
            //     candRR.changeTo(s.SV[i], i);
            // }
            int minPQ2 = std::min(s.p_prime, s.q_prime);
            double zInSubSpace = s.ZCount;

            edges.clear();
            probabilitiesForFirstE.clear();
            probForFirstE.clear();
            aliasForFirstE.clear();

            if (minPQ2 == 1) {
                for (int v = 0; v < s.SVSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ2][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < s.SUSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        // printf("dpU[minPQ2][mapUtoV[i]]: %f\n", dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                        probabilitiesForFirstE.push_back(dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }
            double sum = 0;
            for (double d : probabilitiesForFirstE) {
                sum += d;
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

            for (int i = 0; i < batch.second; i++) {
                t++;
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                int sampledIndexForFirstE = getIndexAlias(gen, probForFirstE, aliasForFirstE);
                edge = edges[sampledIndexForFirstE];

                if (minPQ2 == 1) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                for (int i = 1; i < minPQ2; i++) {
                    // TODO update to use binary search to find the vertex that start the
                    // out neighborhood select u
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ2 - i][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ2 - i + 1][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }

                    // select v
                    if (i != minPQ2 - 1) {
                        temp = 0.0;
                        point = distribution(eng);
                        bool f = false;
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ2 - i][mapUtoV[j]];
                                if (temp + 1e-8 >= point * dpV[minPQ2 - i][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    f = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                int vCount = q - p + 1;
                if (minPQ2 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ2 == 1) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // for (size_t i = 0; i < selectedL.size(); i++) {
                //     printf("%d ", selectedL[i]);
                // }
                // printf("\n");
                assert(selectedR.size() == s.q_prime);
                assert(selectedL.size() == minPQ2);
                // check whether the sample is a biclique
                bool connected = true;

                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candLVec[selectedL[i]], candRVec[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
                successSamples++;

                if (successSamples >= (int)gamma) {
                    gammaReached = true;
                    break;
                }
            }

            if (gammaReached) break;
        }
    }

    printf("pqBicluqeContainZ: %d\n", pqBicluqeContainZ);
    printf("t: %d\n", t);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("successSamples: %d\n", successSamples);
    printf("density of biclique: %f\n", pqBicluqeContainZ / (double)t);
    printf("biclique count: %f\n", pqBicluqeContainZ / (double)t * totalZinShadow);
    printf("time: %f\n", (double)(std::chrono::high_resolution_clock::now() - start).count() / 1e9);
}

// shadow builder with decresinng order

void accuracy::shadowBuilderZStar3(int p, int q, double epsilon) {
    std::signal(SIGSEGV, segfault_handler);
    if (p > q) {
        swap(p, q);
    }
    printf(" Working on p q e: %d - %d - %.2f \n", p, q, epsilon);
    auto start = chrono::high_resolution_clock::now();
    minPQ = min(p, q);

    // get the total structure count in the input graph
    // e1 hols u's neighbors

    std::vector<uint32_t> candLtemp(g->maxDv);
    std::vector<uint32_t> outPosition(g->maxDu);
    vector<int> twoHopCount(g->n1);
    candL.resize(g->n1);
    candR.resize(g->n2);
    std::vector<uint32_t> candLVec, candRVec;
    double totalZinShadow = 0;
    double totalZinShadowTest = 0;

    uint32_t maxE = std::min((1u << 21), g->maxDu * g->maxDv);
    // for Shadow
    uint32_t cnt_PQ = 1;

    priority_queue<Subspace, std::vector<Subspace>, CompareSubspace> shadow;
    vector<Subspace> subspaceTest;
    uint32_t totalS = 0;

    std::vector<uint32_t> pV(g->maxDu);  // Adjusted size to accommodate pV[v + 1]
    std::vector<uint32_t> pU(g->maxDv);
    std::vector<uint32_t> e1(maxE);  // Assuming maximum possible size
    std::vector<uint32_t> e2(
        maxE);  // Assuming maximum possible sizmaxPossibleLengthe
    std::vector<double> ddp(maxE, 0.0);
    std::vector<std::vector<double>> dpU(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<std::vector<double>> dpV(minPQ + 1,
                                         std::vector<double>(maxE, 0.0));
    std::vector<uint32_t> mapUtoV(maxE);
    std::vector<uint32_t> mapVtoU(maxE);
    const double delta = 0.05;
    double gamma = 4 * (1 + epsilon) * (std::exp(1) - 2) * std::log(2 / delta) /
                   (epsilon * epsilon);

    double tSample = 0;
    double tTotal = 0;
    double auxTotalSamples = 1000;
    double pqBicluqeEstAux = 1;
    int candLSum = 0;
    int candRSum = 0;

    for (uint32_t e = 0; e < maxE; e++) {
        dpU[0][e] = 0;
        dpU[1][e] = 1;
        dpV[0][e] = 1;
    }

    // DP Buidling
    auto buildDP = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candL[j];
            pU[j + 1] = pU[j];
            if (rSize < g->deg1(x)) {
                for (int k = 0; k < rSize; k++) {
                    uint32_t y = candR[k];
                    if (g->connectUV(x, y)) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            } else {
                auto start = g->e1.begin() + g->pU[x];
                auto end = g->e1.begin() + g->pU[x + 1];
                uint32_t i = std::upper_bound(start, end, candR[0]) - g->e1.begin();
                for (; i < g->pU[x + 1]; i++) {
                    int k = candR.idx(g->e1[i]);
                    if (k < rSize) {
                        e1[pU[j + 1]++] = k;
                        pV[k + 1]++;
                    }
                }
            }
            // printf("pU[%d]: %d ", j, pU[j]);
        }
        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu

                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    // printf("p: %d, q: %d, C[%d][%d]: %f\n", p, q, pU[u + 1] - i - 1, q - p + 1, C[pU[u + 1] - i - 1][q - p + 1]);
                    dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }
        return k;
    };
    auto buildDP2Vec = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);

        std::fill(pV.begin(), pV.begin() + rSize + 1, 0);

        for (int j = 0; j < lSize; j++) {
            uint32_t x = candLVec[j];
            pU[j + 1] = pU[j];
            for (int k = 0; k < rSize; k++) {
                uint32_t y = candRVec[k];
                if (g->connectUV(x, y)) {
                    e1[pU[j + 1]++] = k;
                    pV[k + 1]++;
                }
            }
        }

        assert(pU[lSize] < maxE);

        if (pU[lSize] == 0) return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0) return 0;
        assert(pU[lSize] == pV[rSize]);

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                int v = e1[i];
                e2[pV[v]] = u;
                mapUtoV[i] = pV[v];
                mapVtoU[pV[v]] = i;
                pV[v]++;
            }
        }
        for (int v = rSize; v >= 1; v--) {
            pV[v] = pV[v - 1];
        }
        pV[0] = 0;

        if (minPQ == 1) {
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u]; i < pU[u + 1]; i++) {
                    // 2 paths for vu
                    if (p == q) {
                        dpV[1][i] = pU[u + 1] - i - 1;
                    } else {
                        dpV[1][i] = C[pU[u + 1] - i - 1][q - p];
                    }
                }
            }
            return 2;
        }
        for (uint32_t u = 0; u < lSize; u++) {
            for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu
                if (p == q) {
                    dpV[1][i] = pU[u + 1] - i - 1;
                } else {
                    seg_value = dpV[1].size();
                    pU_seg = i;
                    pU1_seg = maxE;
                    if (pU[u + 1] - i - 1 == 0) {
                        dpV[1][i] = 0;
                    } else {
                        double y = C[pU[u + 1] - i - 1][q - p + 1];
                        dpV[1][i] = y;
                    }
                    // dpV[1][i] = C[pU[u + 1] - i - 1][q - p + 1];
                }
            }
        }

        int minLR = std::min(lSize, rSize);
        int k = 2;

        // if (minLR < minPQ) return 0;

        for (; k <= minPQ && k <= minLR; k++) {
            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int v = 0; v < rSize; v++) {
                for (int i = pV[v] + 1; i < pV[v + 1]; i++) {
                    if (dpV[k - 1][mapVtoU[i]] < 0.5)
                        continue;
                    ddp[pV[v]] += dpV[k - 1][mapVtoU[i]];
                    ddp[i] -= dpV[k - 1][mapVtoU[i]];
                }
            }
            for (int v = 0; v < rSize; v++) {
                dpU[k][pV[v]] = ddp[pV[v]] < 0 ? 0 : ddp[pV[v]];
                for (int e = pV[v] + 1; e < pV[v + 1]; e++) {
                    dpU[k][e] = dpU[k][e - 1] + ddp[e];
                }
            }

            std::fill(ddp.begin(), ddp.begin() + pU[lSize], 0);
            for (int u = 0; u < lSize; u++) {
                for (int i = pU[u] + 1; i < pU[u + 1]; i++) {
                    if (dpU[k][mapUtoV[i]] < 0.5)
                        continue;
                    ddp[pU[u]] += dpU[k][mapUtoV[i]];
                    ddp[i] -= dpU[k][mapUtoV[i]];
                }
            }
            for (int u = 0; u < lSize; u++) {
                dpV[k][pU[u]] = ddp[pU[u]] < 0 ? 0 : ddp[pU[u]];
                for (int e = pU[u] + 1; e < pU[u + 1]; e++) {
                    dpV[k][e] = dpV[k][e - 1] + ddp[e];
                    // /printf("dpV[%d][%d]: %f\n", k+1, e, dpV[k+1][e]);
                }
            }
        }

        return std::min(k, minLR);
    };

    for (uint32_t u = 0; u < g->n1; u++) {
        if (g->deg1(u) < q) {
            continue;
        }

        uint32_t candRSize = 0;

        // for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
        //     uint32_t v = g->e1[i];
        //     candR.changeTo(v, candRSize);
        //     candRSize++;
        // }

        candLtemp.clear();

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            candR.changeTo(v, candRSize);
            candRSize++;
            for (uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if (w == u) {
                    // keepig the index of u so it can later get vertices larger than u
                    outPosition[i - g->pU[u]] = j;
                    break;
                }
                if (twoHopCount[w] == 0) {
                    candLtemp.push_back(w);
                }
                twoHopCount[w]++;
            }
        }

        int countCandLtemp = 0;
        for (auto w : candLtemp) {
            if (twoHopCount[w] > 1) {
                candL.changeTo(w, countCandLtemp);
                countCandLtemp++;
            }
            twoHopCount[w] = 0;
        }
        // /if no vertices found with two hop neighbors
        if (countCandLtemp == 0) {
            continue;
        }

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            assert(v == candR[0]);
            // if (v != candR[0]) {
            //     printf("v: %d\n", v);
            //     for (int i = 0; i < candRSize; i++) {
            //         printf("candR[%d]: %d\n", i, candR[i]);
            //     }
            // }
            int candLSize = 0;
            // int outPosE = std::find(g->pV[v], g->pV[v + 1], u);
            // auto it = std::find(std::next(g->e2.begin(), g->pV[v]), std::next(g->e2.begin(), g->pV[v + 1]), u);
            // int outPosE = std::distance(g->e2.begin(), it);
            //  /assert(outPosE == outPosition[i - g->pU[u]]);
            for (uint32_t j = outPosition[i - g->pU[u]] + 1; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                // check whether the w is inside the candLTemp
                if (candL.idx(w) < countCandLtemp) candL.changeTo(w, candLSize++);
            }

            candR.changeTo(v, --candRSize);

            // if (candLSize < p - 1)
            //     continue;
            // if (candRSize < q - 1)
            //     break;
            for (int i = 1; i < candRSize; i++) {
                candR.changeToByPos(i, i - 1);
            }
            if (candRSize < q - 1)
                break;
            if (candLSize < p - 1)
                continue;

            // adding to shadow
            // calculate mu here
            double mu = 0;
            // printf("candLSize: %d, candRSize: %d, p: %d, q: %d\n", candLSize,
            // candRSize, p, q); sampling inside shadow to calcule mu
            int maxPossibleLength = buildDP(candLSize, candRSize, p - 1, q - 1);
            if (maxPossibleLength < minPQ - 1) {
                continue;
            }

            double zInSubSpace = 0;
            double zInSubSpaceTest = 0;

            if (minPQ - 1 == 1) {
                for (uint32_t i = 0; i < pV[candRSize]; i++) {
                    zInSubSpace += dpV[minPQ - 1][i];
                }
            } else {
                for (uint32_t i = 0; i < pU[candLSize]; i++) {
                    zInSubSpace += dpU[minPQ - 1][i];
                }
            }
            // if (candL[0] == 555 && candR[0] == 1773) {
            //     printf("success\n");
            //     printf("zInSubSpace: %f\n", zInSubSpace);
            //     printf("candLSize: %d, candRSize: %d\n", candLSize, candRSize);
            //     for (int i = 0; i < candLSize; i++) {
            //         printf("candL[%d]: %d\n", i, candL[i]);
            //     }
            //     for (int i = 0; i < candRSize; i++) {
            //         printf("candR[%d]: %d\n", i, candR[i]);
            //     }
            //     for (int i = 0; i <= candLSize; i++) {
            //         printf("pU[%d]: %d\n", i, pU[i]);
            //     }
            //     for (int i = 0; i <= candRSize; i++) {
            //         printf("pV[%d]: %d\n", i, pV[i]);
            //     }
            //     for (int i = 0; i < pU[candLSize]; i++) {
            //         printf("e1[%d]: %d\n", i, e1[i]);
            //     }
            //     for (int i = 0; i < pV[candRSize]; i++) {
            //         printf("e2[%d]: %d\n", i, e2[i]);
            //     }
            //     exit(0);
            // }
            totalZinShadow += zInSubSpace;

            if (zInSubSpace < 1)
                continue;
            double pqBicluqeContainZ = 0;
            double auxSampleSize = 10;
            int auxSampleSize_temp = auxSampleSize;
            double t_sample = 0;

            auto auxSamplingStart = chrono::high_resolution_clock::now();
            std::pair<uint32_t, uint32_t> edge;
            std::vector<std::pair<uint32_t, uint32_t>> edges;
            vector<double> probabilities;
            vector<double> probabilitiesForFirstE;
            vector<double> Prob, probForFirstE;
            vector<uint32_t> Alias, aliasForFirstE;

            if (minPQ == 2) {
                for (int v = 0; v < candRSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ - 1][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < candLSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpU[minPQ - 1][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

            while (auxSampleSize_temp--) {
                // selecting the first edge
                std::random_device rd;

                std::mt19937 eng(rd());  // Engine (Mersenne Twister)
                std::uniform_real_distribution<> distribution(
                    0.0, 1.0);
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                // if (minPQ - 1 == 1) {
                //     for (int v = 0; v < candRSize; v++) {
                //         for (int i = pV[v]; i < pV[v + 1]; i++) {
                //             temp += dpV[minPQ - 1][mapVtoU[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(v);
                //                 selectedL.push_back(e2[i]);
                //                 preU = e2[i];
                //                 preV = v;
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // } else {
                //     double temp = 0.0;
                //     for (int u = 0; u < candLSize; u++) {
                //         for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                //             temp += dpU[minPQ - 1][mapUtoV[i]];
                //             if (temp + 1e-8 >= point * zInSubSpace) {
                //                 selectedR.push_back(e1[i]);
                //                 selectedL.push_back(u);
                //                 preU = u;
                //                 preV = e1[i];
                //                 preE = i;
                //                 selectedLEmpty = false;
                //                 break;
                //             }
                //         }
                //         if (!selectedLEmpty) {
                //             break;
                //         }
                //     }
                // }

                int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                edge = edges[sampledIndexForFirstE];
                if (minPQ == 2) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                // printf("minPQ: %d\n", minPQ);
                // printf("---------------------\n");
                // printf("preV: %d, preU: %d\n", preV, preU);
                for (int i = 1; i < minPQ - 1; i++) {
                    // out neighborhood select u
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ - i - 1][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ - i][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }
                    // select v
                    if (i != minPQ - 2) {
                        // printf("dpV[%d][%d]: %f\n", minPQ - i - 1, mapVtoU[preE], dpV[minPQ - i - 1][mapVtoU[preE]]);
                        // printf("pU[preU]: %d, pU[preU + 1]: %d\n", pU[preU], pU[preU + 1]);
                        // printf("preV: %d\n", preV);
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ - i - 1][mapUtoV[j]];
                                // printf("temp: %f\n", temp);
                                // printf("point * dpV[minPQ - i - 1][mapVtoU[preE]]: %f\n", point * dpV[minPQ - i - 1][mapVtoU[preE]]);
                                if (temp + 1e-8 >= point * dpV[minPQ - i - 1][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                    }
                }
                if (selectedL.size() != minPQ - 1) {
                    printf("selectedL.size(): %d, minPQ - 1: %d\n", selectedL.size(), minPQ - 1);
                    continue;
                }
                assert(selectedL.size() == minPQ - 1);
                int vCount = q - p + 1;
                if (minPQ - 1 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ == 2) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                // assert(pU[preU + 1] - index >= vCount);
                // assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // if (selectedR.size() != q - 1) {
                //     printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                //     printf("selectedL.size(): %d, p - 1: %d\n", selectedL.size(), p - 1);
                //     exit(0);
                // }
                if (selectedR.size() != q - 1) {
                    printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                    continue;
                }
                assert(selectedR.size() == q - 1);
                assert(selectedL.size() == minPQ - 1);

                // check whether the sampled zstar is a biclique
                bool connected = true;
                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candL[selectedL[i]], candR[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
            }

            // calculating t_sample
            auto auxSamplingEnd = chrono::high_resolution_clock::now();
            double auxSamplingDuration = chrono::duration_cast<chrono::microseconds>(
                                             auxSamplingEnd - auxSamplingStart)
                                             .count();
            tTotal += auxSamplingDuration;
            auxTotalSamples += auxSampleSize;
            // printf("pqBicluqeContainZ: %f\n", pqBicluqeContainZ   );
            pqBicluqeEstAux += zInSubSpace * (pqBicluqeContainZ / auxSampleSize);

            // subspaceTest.push_back(Subspace(p - 1, q - 1, candL, candR, candLSize,
            // candRSize, pqBicluqeContainZ / auxSampleSize, zInSubSpace *
            // (pqBicluqeContainZ / auxSampleSize), zInSubSpace));
            std::vector<uint32_t> candLShadow;
            std::vector<uint32_t> candRShadow;
            for (int i = 0; i < candLSize; i++) {
                candLShadow.push_back(candL[i]);
            }
            for (int i = 0; i < candRSize; i++) {
                candRShadow.push_back(candR[i]);
            }

            if (zInSubSpace > 0.5) {
                shadow.push(Subspace(p - 1, q - 1, candLShadow, candRShadow, candLSize, candRSize,
                                     pqBicluqeContainZ / auxSampleSize,
                                     zInSubSpace * (pqBicluqeContainZ / auxSampleSize),
                                     zInSubSpace));
            }

            t_sample += auxSamplingDuration;

            // if (maxPossibleLength >= minPQ) {
            // 	if (pU[minPQ - 1] > 0.3) {
            // 		for (int i = 0; i < candLSize; i++) {
            // 			totalZinShadow += dpU[minPQ - 1][i];
            // 		}
            // 	}
            // }
        }
    }

    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage I is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////
    printf("totalZinShadow: %f\n", totalZinShadow);
    auto end = chrono::high_resolution_clock::now();
    double shadowDuration =
        chrono::duration_cast<chrono::microseconds>(end - start).count();
    std::random_device rd;
    std::mt19937 eng(rd());  // Engine (Mersenne Twister)
    std::uniform_real_distribution<> distribution(
        0.0, 1.0);  // Distribution for doubles between 0.0 and 1.0 (inclusive of
    // 0.0 and exclusive of 1.0)
    vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
    printf("gamma: %f\n", gamma);
    printf("tTotal: %f\n", tTotal);
    printf("auxTotalSamples: %f\n", auxTotalSamples);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("Biclique Density: %f\n", pqBicluqeEstAux / totalZinShadow);
    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    fflush(stdout);
    printf("shadow size before the second refinement : %d \n", shadow.size());
    printf(" gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux: %f\n", gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    int count = 0;
    printf("shadowDuration: %f\n", shadowDuration);
    printf("shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux : %d\n", shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux);
    // /shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux
    fflush(stdout);
    if (shadow.size() == 0) {
        printf("No biclique found fopr p: %d, q: %d\n", p, q);
        return;
    }

    while (shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux) {
        Subspace min_space = shadow.top();
        if (min_space.mu > 0.5) break;
        shadow.pop();
        if (min_space.p_prime == 1 || min_space.q_prime == 1) {
            min_space.mu = 1;
            shadow.push(min_space);
        }

        uint32_t llSize = min_space.SUSize;
        uint32_t rrSize = min_space.SVSize;
        std::vector<uint32_t> candLL, candRR;
        candLL = min_space.SU;
        candRR = min_space.SV;
        int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
        vector<uint32_t> pUU(g->maxDv + 1);
        vector<uint32_t> pVV(g->maxDu + 1);
        vector<uint32_t> e11(g->m);
        vector<uint32_t> e22(g->m);
        // vector<uint32_t> candLtemp2(g->maxDv * g->maxDv);
        // vector<uint32_t> outPosition2(g->maxDu * g->maxDu);
        std::vector<uint32_t> minSU = min_space.SU;
        std::vector<uint32_t> minSV = min_space.SV;
        std::unordered_map<uint32_t, uint32_t> minSUMap;
        std::unordered_map<uint32_t, uint32_t> minSVMap;
        // vector<uint32_t> twoHopCount2(g->maxDu * g->maxDu);

        pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
        totalZinShadow -= min_space.ZCount;

        for (int i = 0; i < llSize; i++) {
            minSUMap[minSU[i]] = i;
        }
        for (int i = 0; i < rrSize; i++) {
            minSVMap[minSV[i]] = i;
        }

        /////////////////////////////
        /* ***subgraph creating*** */
        ///////////////////////////
        int edgeIndex = 0;
        for (int j = 0; j < llSize; j++) {
            uint32_t x = candLL[j];
            pUU[j + 1] = pUU[j];
            for (int k = 0; k < rrSize; k++) {
                uint32_t y = candRR[k];
                if (g->connectUV(x, y)) {
                    e11[edgeIndex++] = y;
                    pUU[j + 1]++;
                    pVV[k + 1]++;
                }
            }
        }

        // Resize e11 to the actual number of edges
        // Convert pVV to prefix sum array
        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }
        assert(pVV[rrSize] == pUU[llSize]);
        // for (int i = 0; i <= llSize; i++) {
        //     printf("pUU[%d]: %d\n", i, pUU[i]);
        // }
        // for (int i = 0; i <= rrSize; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        // Filling e22
        // edgeIndex = 0;

        // Filling e22
        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                int index = minSVMap[v];
                e22[pVV[index]] = candLL[u];
                pVV[index]++;
            }
        }

        // Reset pVV to prefix sum array
        for (int i = rrSize + 1; i >= 1; i--) {
            pVV[i] = pVV[i - 1];
        }
        pVV[0] = 0;
        // if (pUU[llSize] != pVV[rrSize]) {
        //     printf("minSpace.p_prime: %d, minSpace.q_prime: %d\n", min_space.p_prime, min_space.q_prime);
        //     printf("--------------------\n");

        //     for (int i = 0; i <= llSize; i++) {
        //         printf("pUU[%d]: %d\n", i, pUU[i]);
        //     }
        //     for (int i = 0; i <= rrSize; i++) {
        //         printf("pVV[%d]: %d\n", i, pVV[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < pUU[llSize]; i++) {
        //         printf("e11[%d]: %d\n", i, e11[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < pVV[rrSize]; i++) {
        //         printf("e22[%d]: %d\n", i, e22[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < llSize; i++) {
        //         printf("candLL[%d]: %d\n", i, candLL[i]);
        //     }
        //     printf("--------------------\n");
        //     for (int i = 0; i < rrSize; i++) {
        //         printf("candRR[%d]: %d\n", i, candRR[i]);
        //     }
        //     exit(0);
        // }
        assert(pVV[rrSize] == pUU[llSize]);

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////
        // printf("min_space.subgraph: %d\n", min_space.subgraph);
        // printf("min_space.p_prime: %d\n", min_space.p_prime);
        // for (int i = 0; i < llSize; i++) {
        //     printf("candLL[%d]: %d\n", i, candLL[i]);
        // }
        // for (int i = 0; i < rrSize; i++) {
        //     printf("candRR[%d]: %d\n", i, candRR[i]);
        // }
        // for (int i = 0; i <= llSize; i++) {
        //     printf("pUU[%d]: %d\n", i, pUU[i]);
        // }
        // for (int i = 0; i <= rrSize; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        // for (int i = 0; i <= pUU[llSize]; i++) {
        //     printf("e11[%d]: %d\n", i, e11[i]);
        // }
        // for (int i = 0; i <= pVV[rrSize]; i++) {
        //     printf("e22[%d]: %d\n", i, e22[i]);
        // }
        candLL.clear();
        candRR.clear();

        uint32_t totalZinThisSub = 0;
        for (uint32_t uIndex = 0; uIndex < llSize; uIndex++) {
            candRR.clear();
            uint32_t u = minSU[uIndex];
            uint32_t candRSize2 = 0;

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                uint32_t v = e11[i];
                candRR.push_back(v);
                candRSize2++;
            }
            // for (int i = 0; i < candRSize2; i++) {
            //     printf(" main loop candRR[%d]: %d\n", i, candRR[i]);
            // }
            // candLtemp2.clear();

            // for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
            //     uint32_t v = e11[i];
            //     uint32_t vIndex = minSVMap[v];
            //     for (uint32_t j = pVV[vIndex + 1] - 1; j >= pVV[vIndex]; j--) {
            //         uint32_t w = e22[j];
            //         if (w == u) {
            //             // keepig the index of u so it can later get vertices larger than u
            //             outPosition2[i - pUU[uIndex]] = j;
            //             break;
            //         }
            //         // if (twoHopCount2[w] == 0) {
            //         //     candLtemp2.push_back(w);
            //         // }
            //         // twoHopCount2[w]++;
            //     }
            // }

            // int countCandLtemp2 = 0;
            // for (auto w : candLtemp2) {
            //     if (twoHopCount2[w] > 0) {
            //         candLL.push_back(w);
            //         countCandLtemp2++;
            //     }
            //     twoHopCount2[w] = 0;
            // }

            // // if no vertices found with two hop neighbors

            // if (countCandLtemp2 == 0) {
            //     continue;
            // }

            for (uint32_t i = pUU[uIndex]; i < pUU[uIndex + 1]; i++) {
                candLL.clear();
                uint32_t v = e11[i];
                assert(v == candRR[0]);
                int candLSize2 = 0;

                // uint32_t j = outPosition2[i - pUU[uIndex]] + 1;
                uint32_t j = std::find(e22.begin() + pVV[minSVMap[v]], e22.begin() + pVV[minSVMap[v] + 1], u) - e22.begin() + 1;
                uint32_t vIndex = minSVMap[v];

                for (; j < pVV[vIndex + 1]; j++) {
                    uint32_t w = e22[j];
                    // check whether the w is inside the candL
                    candLL.push_back(w);
                    candLSize2++;
                }
                // fflush(stdout);
                // if (candRR[0] != v) {
                //     printf("candR is not equal to v\n");
                //     printf("v: %d\n", v);
                //     for (int i = 0; i < candRSize2; i++) {
                //         printf("candRR[%d]: %d\n", i, candRR[i]);
                //     }
                //     fflush(stdout);
                //     exit(0);
                // }

                //  candRR.changeTo(v, --candRSize2);
                // printf("candR befdore: %d\n", candRSize2);
                // for (int i = 0; i < candRSize2; i++) {
                //     printf("candRR[%d]: %d\n", i, candRR[i]);
                // }
                candRSize2--;
                candRR.erase(candRR.begin());
                // printf("candR after: %d\n", candRSize2);
                // for (int i = 0; i < candRSize2; i++) {
                //     printf("candRR[%d]: %d\n", i, candRR[i]);
                // }
                if (candRSize2 < min_space.q_prime - 1) {
                    break;
                }
                if (candLSize2 < min_space.p_prime - 1) {
                    continue;
                }

                // for (int i = 0; i < candRSize2; i++) {
                //     candRR.changeToByPos(i, i - 1);
                // }

                candLVec = candLL;
                candRVec = candRR;

                int maxPossibleLength = buildDP2Vec(candLSize2, candRSize2, min_space.p_prime - 1, min_space.q_prime - 1);

                if (maxPossibleLength < minPQ2 - 1) {
                    continue;
                }
                double zInMinSubSpace = 0;
                if (minPQ2 == 2) {
                    for (uint32_t i = 0; i < pV[candRSize2]; i++) {
                        zInMinSubSpace += dpV[minPQ2 - 1][i];
                    }

                }

                else {
                    for (uint32_t i = 0; i < pU[candLSize2]; i++) {
                        zInMinSubSpace += dpU[minPQ2 - 1][i];
                    }
                }

                totalZinThisSub += zInMinSubSpace;
                totalZinShadow += zInMinSubSpace;

                if (zInMinSubSpace < 1)
                    continue;

                // sampling
                uint32_t auxSampleSize = 100;
                uint32_t auxSampleSizeTemp = auxSampleSize;
                uint32_t pqBicluqeContainZ = 0;
                std::pair<uint32_t, uint32_t> edge;
                std::vector<std::pair<uint32_t, uint32_t>> edges;
                vector<double> probabilities;
                vector<double> probabilitiesForFirstE;
                vector<double> Prob, probForFirstE;
                vector<uint32_t> Alias, aliasForFirstE;

                if (minPQ2 == 2) {
                    for (int v = 0; v < candRSize2; v++) {
                        for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                            edge.first = i;
                            edge.second = v;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpV[minPQ2 - 1][mapVtoU[i]] / zInMinSubSpace);
                        }
                    }
                } else {
                    for (int u = 0; u < candLSize2; u++) {
                        for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                            edge.first = u;
                            edge.second = i;
                            edges.push_back(edge);
                            probabilitiesForFirstE.push_back(dpU[minPQ2 - 1][mapUtoV[i]] / zInMinSubSpace);
                        }
                    }
                }

                initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

                while (auxSampleSizeTemp--) {
                    // selecting the first edge
                    selectedL.clear();
                    selectedR.clear();
                    double point = distribution(eng);
                    double temp = 0.0;
                    uint32_t preU = 0, preV = 0, preE = 0;
                    bool selectedLEmpty = true;

                    // //pasted
                    // for (int u = 0; u < candLSize2; u++) {
                    // 	for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                    // 		temp += dpU[minPQ2 - 1][mapUtoV[i]];
                    // 		if (temp + 1e-8 >= point * zInMinSubSpace) {
                    // 			selectedR.push_back(e1[i]);
                    // 			selectedL.push_back(u);
                    // 			preU = u;
                    // 			preV = e1[i];
                    // 			preE = i;
                    // 			selectedLEmpty = false;
                    // 			break;
                    // 		}
                    // 	}
                    // 	if (!selectedLEmpty) {
                    // 		break;
                    // 	}
                    // }
                    int sampledIndexForFirstE = getIndexAlias(eng, probForFirstE, aliasForFirstE);

                    edge = edges[sampledIndexForFirstE];
                    if (minPQ2 == 2) {
                        selectedL.push_back(e2[edge.first]);
                        selectedR.push_back(edge.second);
                        preE = edge.first;
                        preU = e2[edge.first];
                        preV = edge.second;
                    } else {
                        selectedR.push_back(e1[edge.second]);
                        selectedL.push_back(edge.first);
                        preE = edge.second;
                        preU = edge.first;
                        preV = e1[edge.second];
                    }

                    assert(selectedR.size() > 0);
                    assert(selectedL.size() > 0);
                    // selecting the rest of the edge
                    for (int i = 1; i < minPQ2 - 1; i++) {
                        // out neighborhood select u
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                            if (e2[j] > preU) {
                                temp += dpV[minPQ2 - i - 1][mapVtoU[j]];
                                if (temp + 1e-8 >= point * dpU[minPQ2 - i][mapUtoV[preE]]) {
                                    uint32_t u = e2[j];
                                    selectedL.push_back(u);
                                    preU = u;
                                    preE = j;
                                    break;
                                }
                            }
                        }
                        // select v
                        if (i != minPQ2 - 2) {
                            temp = 0.0;
                            point = distribution(eng);
                            for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                                if (e1[j] > preV) {
                                    temp += dpU[minPQ2 - i - 1][mapUtoV[j]];
                                    if (temp + 1e-8 >= point * dpV[minPQ2 - i - 1][mapVtoU[preE]]) {
                                        uint32_t v = e1[j];
                                        selectedR.push_back(v);
                                        preV = v;
                                        preE = j;
                                        break;
                                    }
                                }
                            }
                        }
                    }

                    int vCount = q - p + 1;
                    if (minPQ2 == 2) {
                        vCount = q - p;
                    }

                    uint32_t index = 0;
                    if (minPQ == 2) {
                        index = pU[preU];
                    } else {
                        auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                        index = std::distance(e1.begin(), it);
                    }

                    assert(index >= pU[preU] && index < pU[preU + 1]);

                    std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                    for (int i = 0; i < selectedRStar.size(); i++) {
                        selectedR.push_back(e1[selectedRStar[i]]);
                    }
                    if (selectedR.size() != min_space.q_prime - 1 || selectedL.size() != minPQ2 - 1) {
                        printf("issue is sample zsatar sizes \n");
                        printf("selectedR.size(): %d, q - 1: %d\n", selectedR.size(), q - 1);
                        continue;
                    }
                    assert(selectedR.size() == min_space.q_prime - 1);
                    assert(selectedL.size() == minPQ2 - 1);
                    // check whether the sampled z is a biclique
                    bool connected = true;

                    for (int i = 0; i < selectedL.size(); i++) {
                        for (int j = 0; j < selectedR.size(); j++) {
                            if (i == j)
                                continue;
                            if (i > 0 && i == j + 1)
                                continue;
                            if (!g->connectUV(candLVec[selectedL[i]], candRVec[selectedR[j]])) {
                                connected = false;
                                break;
                            }
                        }
                        if (!connected)
                            break;
                    }
                    if (!connected)
                        continue;
                    pqBicluqeContainZ++;
                }
                Subspace s = Subspace(min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR, candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize, zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize), zInMinSubSpace);
                s.subgraph = s.subgraph + 1;
                shadow.push(s);
                // shadow.push(Subspace(
                //     min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR,
                //     candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize,
                //     zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize),
                //     zInMinSubSpace));

                pqBicluqeEstAux +=
                    zInMinSubSpace * pqBicluqeContainZ / (double)auxSampleSize;
            }
        }

        end = std::chrono::high_resolution_clock::now();
        shadowDuration =
            std::chrono::duration_cast<std::chrono::microseconds>(end - start)
                .count();

        candLL.clear();
        candRR.clear();
        e11.clear();
        e22.clear();
        pVV.clear();
        // outPosition2.clear();
        // twoHopCount2.clear();
    }

    printf("end the second phrase\n");
    printf("shadow size: %d\n", shadow.size());

    printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
    printf("totalZinShadow: %f\n", totalZinShadow);
    double roughPQEstDens = pqBicluqeEstAux / totalZinShadow;
    vector<Subspace> subs;
    uint32_t batchSize = gamma / roughPQEstDens;
    std::pair<uint32_t, uint32_t> edge;
    std::vector<std::pair<uint32_t, uint32_t>> edges;
    printf("roughPQEstDens: %f\n", roughPQEstDens);
    vector<double> probabilities;
    vector<double> probabilitiesForFirstE;
    vector<Subspace> shadowVec;
    vector<double> Prob, probForFirstE;
    vector<uint32_t> Alias, aliasForFirstE;
    std::mt19937 gen(rd());
    uint32_t pqBicluqeContainZ = 0;

    printf("shadow.size(): %d\n", shadow.size());
    if (shadow.size() == 0) {
        printf("No biclique found fopr p: %d, q: %d\n", p, q);
        return;
    }
    /////////////////////////////////////////////
    /////////////////////////////////////////////
    /* ***Shadow Refinement stage II is done*** */
    ////////////////////////////////////////////
    ////////////////////////////////////////////

    while (!shadow.empty()) {
        shadowVec.push_back(shadow.top());
        shadow.pop();
    }

    printf("totalZinShadow: %f\n", totalZinShadow);

    for (const auto& subspace : shadowVec) {
        probabilities.push_back(static_cast<double>(subspace.ZCount) / totalZinShadow);
    }

    initializeAliasMethod(probabilities, Prob, Alias);
    int sampledIndex = getIndexAlias(gen, Prob, Alias);
    std::unordered_map<int, int> subSpaceAndC;
    printf("batchSize: %d\n", batchSize);
    int successSamples = 0;
    int t = 0;
    bool gammaReached = false;
    double totalZFinal = 0;

    while (successSamples < int(gamma)) {
        subSpaceAndC.clear();

        for (int i = 0; i < batchSize; i++) {
            int sampledIndex = getIndexAlias(gen, Prob, Alias);
            if (subSpaceAndC.count(sampledIndex) == 0) {
                subSpaceAndC[sampledIndex] = 1;
            } else {
                subSpaceAndC[sampledIndex] = subSpaceAndC[sampledIndex] + 1;
            }
        }

        for (const auto& batch : subSpaceAndC) {
            Subspace s = shadowVec[batch.first];
            totalZFinal += s.ZCount;
            candLVec.clear();
            candRVec.clear();
            // candL.resize(g->n1 + 1);
            // candR.resize(g->n2 + 1);
            candLVec = s.SU;
            candRVec = s.SV;
            // printf("s.SU.size(): %d\n", s.SU.size());
            // printf("s.SUSize: %d\n", s.SUSize);
            // for (int i = 0; i < s.SUSize; i++) {
            //     printf("s.SU[i]: %d\n", s.SU[i]);
            // }
            // continue;
            // exit(0);
            // for (int i = 0; i < s.SUSize; i++) {
            //     candL.changeTo(s.SU[i], i);
            // }
            // for (int i = 0; i < s.SVSize; i++) {
            //     candR.changeTo(s.SV[i], i);
            // }

            int maxLength = buildDP2Vec(s.SUSize, s.SVSize, s.p_prime, s.q_prime);
            uint32_t llSize = s.SUSize;
            uint32_t rrSize = s.SVSize;
            // LinearSet candLL, candRR;
            // candLL.resize(g->n1 + 1);
            // candRR.resize(g->n2 + 1);
            // candLL = s.SU;
            // candRR = s.SV;
            // for (int i = 0; i < s.SUSize; i++) {
            //     candLL.changeTo(s.SU[i], i);
            // }
            // for (int i = 0; i < s.SVSize; i++) {
            //     candRR.changeTo(s.SV[i], i);
            // }
            int minPQ2 = std::min(s.p_prime, s.q_prime);
            double zInSubSpace = s.ZCount;

            edges.clear();
            probabilitiesForFirstE.clear();
            probForFirstE.clear();
            aliasForFirstE.clear();

            if (minPQ2 == 1) {
                for (int v = 0; v < s.SVSize; v++) {
                    for (uint32_t i = pV[v]; i < pV[v + 1]; i++) {
                        edge.first = i;
                        edge.second = v;
                        edges.push_back(edge);
                        probabilitiesForFirstE.push_back(dpV[minPQ2][mapVtoU[i]] / zInSubSpace);
                    }
                }
            } else {
                for (int u = 0; u < s.SUSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        edge.first = u;
                        edge.second = i;
                        edges.push_back(edge);
                        // printf("dpU[minPQ2][mapUtoV[i]]: %f\n", dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                        probabilitiesForFirstE.push_back(dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                    }
                }
            }
            double sum = 0;
            for (double d : probabilitiesForFirstE) {
                sum += d;
            }

            initializeAliasMethod(probabilitiesForFirstE, probForFirstE, aliasForFirstE);

            for (int i = 0; i < batch.second; i++) {
                t++;
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                int sampledIndexForFirstE = getIndexAlias(gen, probForFirstE, aliasForFirstE);
                edge = edges[sampledIndexForFirstE];

                if (minPQ2 == 1) {
                    selectedL.push_back(e2[edge.first]);
                    selectedR.push_back(edge.second);
                    preE = edge.first;
                    preU = e2[edge.first];
                    preV = edge.second;
                } else {
                    selectedR.push_back(e1[edge.second]);
                    selectedL.push_back(edge.first);
                    preE = edge.second;
                    preU = edge.first;
                    preV = e1[edge.second];
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);

                // selecting the rest of the edge
                for (int i = 1; i < minPQ2; i++) {
                    // TODO update to use binary search to find the vertex that start the
                    // out neighborhood select u
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pV[preV]; j < pV[preV + 1]; j++) {
                        if (e2[j] > preU) {
                            temp += dpV[minPQ2 - i][mapVtoU[j]];
                            if (temp + 1e-8 >= point * dpU[minPQ2 - i + 1][mapUtoV[preE]]) {
                                uint32_t u = e2[j];
                                selectedL.push_back(u);
                                preU = u;
                                preE = j;
                                break;
                            }
                        }
                    }

                    // select v
                    if (i != minPQ2 - 1) {
                        temp = 0.0;
                        point = distribution(eng);
                        bool f = false;
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ2 - i][mapUtoV[j]];
                                if (temp + 1e-8 >= point * dpV[minPQ2 - i][mapVtoU[preE]]) {
                                    uint32_t v = e1[j];
                                    selectedR.push_back(v);
                                    preV = v;
                                    preE = j;
                                    f = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                int vCount = q - p + 1;
                if (minPQ2 == 1) {
                    vCount = q - p;
                }

                uint32_t index = 0;
                if (minPQ2 == 1) {
                    index = pU[preU];
                } else {
                    auto it = std::find(e1.begin() + pU[preU], e1.begin() + pU[preU + 1], preV);
                    index = std::distance(e1.begin(), it);
                }

                assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
                }
                // for (size_t i = 0; i < selectedL.size(); i++) {
                //     printf("%d ", selectedL[i]);
                // }
                // printf("\n");
                assert(selectedR.size() == s.q_prime);
                assert(selectedL.size() == minPQ2);
                // check whether the sample is a biclique
                bool connected = true;

                for (int i = 0; i < selectedL.size(); i++) {
                    for (int j = 0; j < selectedR.size(); j++) {
                        if (i == j)
                            continue;
                        if (i > 0 && i == j + 1)
                            continue;
                        if (!g->connectUV(candLVec[selectedL[i]], candRVec[selectedR[j]])) {
                            connected = false;
                            break;
                        }
                    }
                    if (!connected)
                        break;
                }
                if (!connected)
                    continue;
                pqBicluqeContainZ++;
                successSamples++;

                if (successSamples >= (int)gamma) {
                    gammaReached = true;
                    break;
                }
            }

            if (gammaReached) break;
        }
    }

    printf("pqBicluqeContainZ: %d\n", pqBicluqeContainZ);
    printf("t: %d\n", t);
    printf("totalZinShadow: %f\n", totalZinShadow);
    printf("successSamples: %d\n", successSamples);
    printf("density of biclique: %f\n", pqBicluqeContainZ / (double)t);
    printf("biclique count: %f\n", pqBicluqeContainZ / (double)t * totalZinShadow);
    printf("time: %f\n", (double)(std::chrono::high_resolution_clock::now() - start).count() / 1e9);
}

void accuracy::initializeAliasMethod(vector<double>& probabilities, vector<double>& Prob, vector<uint32_t>& Alias) {
    uint16_t n = probabilities.size();
    Prob.resize(n);
    Alias.resize(n);

    std::vector<double> scaledProbabilities(n);
    std::queue<int> Small, Large;

    for (size_t i = 0; i < n; ++i) {
        scaledProbabilities[i] = probabilities[i] * n;
        if (scaledProbabilities[i] < 1.0) {
            Small.push(i);
        } else {
            Large.push(i);
        }
    }

    while (!Small.empty() && !Large.empty()) {
        int l = Small.front();
        Small.pop();
        int g = Large.front();
        Large.pop();

        Prob[l] = scaledProbabilities[l];
        Alias[l] = g;
        scaledProbabilities[g] =
            (scaledProbabilities[g] + scaledProbabilities[l]) - 1.0;

        if (scaledProbabilities[g] < 1.0) {
            Small.push(g);
        } else {
            Large.push(g);
        }
    }

    while (!Large.empty()) {
        int g = Large.front();
        Large.pop();
        Prob[g] = 1.0;
    }

    while (!Small.empty()) {
        int l = Small.front();
        Small.pop();
        Prob[l] = 1.0;
    }
}

int accuracy::getIndexAlias(std::mt19937& gen, vector<double>& Prob,
                            vector<uint32_t>& Alias) {
    std::uniform_int_distribution<> die(0, Prob.size() - 1);
    int i = die(gen);
    std::uniform_real_distribution<> coin(0.0, 1.0);
    if (coin(gen) < Prob[i]) {
        return i;
    } else {
        return Alias[i];
    }
}

std::vector<uint32_t> accuracy::reservoirSample(std::vector<uint32_t>& vec, int n) {
    std::vector<uint32_t> reservoir(n);
    std::random_device rd;
    std::mt19937 gen(rd());

    // Fill the reservoir array with the first n elements of vec
    for (int i = 0; i < n; ++i) {
        reservoir[i] = vec[i];
    }

    // Replace elements with gradually decreasing probability
    for (int i = n; i < vec.size(); ++i) {
        std::uniform_int_distribution<> dist(0, i);
        int j = dist(gen);
        if (j < n) {
            reservoir[j] = vec[i];
        }
    }

    return reservoir;
}

std::vector<uint32_t> accuracy::sampeleRest(int r, int preU, std::vector<uint32_t>& pU, uint32_t outIndex) {
    std::unordered_set<uint32_t> selectedIndices;
    std::vector<uint32_t> result;
    result.reserve(r);
    std::random_device rd;
    std::mt19937 gen(rd());

    if (outIndex + 1 > pU[preU + 1] - 1) {
        printf("outIndex: %d\n", outIndex + 1);
        printf("pU[preU + 1]: %d\n", pU[preU + 1] - 1);
        exit(0);
    }
    std::uniform_int_distribution<> dis(outIndex + 1, pU[preU + 1] - 1);
    uint32_t selectedIndicesSize = 0;
    // /assert(outIndex + 1 < pU[preU + 1] - 1);
    while (selectedIndicesSize < r) {
        uint32_t index = dis(gen);
        if (selectedIndices.find(index) == selectedIndices.end()) {
            selectedIndices.insert(index);
            result.push_back(index);
            selectedIndicesSize++;
        }
    }

    return result;
}
