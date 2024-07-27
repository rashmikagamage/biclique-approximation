#include "accuracy.h"

#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
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

struct Subspace {
    int p_prime, q_prime;
    LinearSet SU, SV;
    uint32_t SUSize, SVSize;
    double mu;
    double countPQMean;
    double ZCount;
    Subspace(int p_prime, int q_prime, LinearSet SU, LinearSet SV,
             uint32_t SUSize, uint32_t SVSize, double mu, double countPQMean,
             double ZCount)
        : p_prime(p_prime), q_prime(q_prime), SU(SU), SV(SV), SUSize(SUSize), SVSize(SVSize), mu(mu), countPQMean(countPQMean), ZCount(ZCount) {}
};

struct CompareSubspace {
    bool operator()(const Subspace& a, const Subspace& b) { return a.mu > b.mu; }
};

void accuracy::shadowBuilder1(int p, int q, double e) {
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

    uint32_t maxE = std::min((1u << 21), g->maxDu * g->maxDv);
    // for Shadow
    uint32_t cnt_PQ = 1;

    // TODO use fibonacci heap
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

        if (pU[lSize] == 0)

            return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }

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

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu
                dpV[1][i] = pU[u + 1] - i - 1;
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

    auto buildDP2 = [&](int lSize, int rSize, int p, int q) -> int {
        int minPQ = min(p, q);
        pV.clear();
        pU.clear();

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

        if (pU[lSize] == 0)

            return 0;

        for (int v = 0; v < rSize; v++) {
            pV[v + 1] += pV[v];
        }
        if (pV[rSize] == 0)
            return 0;
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

        for (int u = 0; u < lSize; u++) {
            for (int i = pU[u]; i < pU[u + 1]; i++) {
                // 2 paths for vu
                dpV[1][i] = pU[u + 1] - i - 1;
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
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < lSize; j++) {
                printf("dpU[%d][%d]: %f\n", i, j, dpU[i][j]);
            }
        }
        return k;
    };

    for (uint32_t u = 0; u < g->n1; u++) {
        if (g->deg1(u) < 2) {
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
            if (maxPossibleLength < minPQ)
                continue;

            double zInSubSpace = 0;
            for (uint32_t i = 0; i < pU[candLSize]; i++) {
                zInSubSpace += dpU[minPQ - 1][i];
            }
            totalZinShadow += zInSubSpace;

            if (zInSubSpace < 1)
                continue;

            double pqBicluqeContainZ = 0;
            double auxSampleSize = 20;
            int auxSampleSize_temp = auxSampleSize;
            double t_sample = 0;

            auto auxSamplingStart = chrono::high_resolution_clock::now();

            while (auxSampleSize_temp--) {
                // selecting the first edge
                std::random_device rd;

                std::mt19937 eng(rd());  // Engine (Mersenne Twister)
                std::uniform_real_distribution<> distribution(
                    0.0, 1.0);  // Distribution for doubles between 0.0 and 1.0
                // (inclusive of 0.0 and exclusive of 1.0)
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                int in = 0;
                for (int u = 0; u < candLSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        temp += dpU[minPQ - 1][mapUtoV[i]];
                        if (temp + 1e-8 >= point * zInSubSpace) {
                            in = i;
                            selectedR.push_back(e1[i]);
                            selectedL.push_back(u);
                            preU = u;
                            preV = e1[i];
                            preE = i;
                            selectedLEmpty = false;
                            break;
                        }
                    }
                    if (!selectedLEmpty) {
                        break;
                    }
                }

                assert(selectedR.size() > 0);
                assert(selectedL.size() > 0);
                // selecting the rest of the edge
                for (int i = 1; i < minPQ - 1; i++) {
                    // TODO update to use binary search to find the vertex that start the
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
                    temp = 0.0;
                    point = distribution(eng);
                    for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                        if (e1[j] > preV) {
                            temp += dpU[minPQ - i - 1][mapUtoV[j]];
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

                assert(selectedR.size() == minPQ - 1);
                assert(selectedL.size() == minPQ - 1);

                // check whether the sampled z is a biclique
                bool connected = true;
                for (int i = 0; i < minPQ - 1; i++) {
                    for (int j = 0; j < minPQ - 1; j++) {
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
                uint32_t noOfVsConnected = 0;
                if (connected && p == q) {
                    noOfVsConnected = 1;
                }

                for (int i = preV + 1; i < candRSize; i++) {
                    uint32_t v = candR[i];
                    bool check = true;
                    for (int i = 0; i < minPQ - 1; i++) {
                        if (!g->connectUV(candL[selectedL[i]], v)) {
                            check = false;
                            break;
                        }
                    }
                    if (check)
                        noOfVsConnected++;
                }
                if (noOfVsConnected >= q - p) {
                    pqBicluqeContainZ += C[noOfVsConnected][q - p];
                }
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

            shadow.push(Subspace(p - 1, q - 1, candL, candR, candLSize, candRSize,
                                 pqBicluqeContainZ / auxSampleSize,
                                 zInSubSpace * (pqBicluqeContainZ / auxSampleSize),
                                 zInSubSpace));

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

    // getting k paths for each edge in G is done

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

    // fixme working here
    LinearSet sut;
    LinearSet svt;
    sut.resize(g->n1 + 1);
    svt.resize(g->n2 + 1);

    // for (int i = 0; i < g->n1 - 1; i++) {
    //     int val = i * 10 + 1;
    //     if (i * 2 < g->n1) {
    //         sut.changeTo(i * 2, i / 2);
    //     }
    // }
    // for (int i = 0; i < g->n2 - 1; i++) {
    //     int val = i * 10 + 1;
    //     if (i * 2 < g->n2) svt.changeTo(i * 2, i / 2);
    // }

    // sut.changeTo(3, 0);
    // sut.changeTo(6, 1);
    // sut.changeTo(9, 2);
    // svt.changeTo(3, 0);
    // svt.changeTo(6, 1);
    // svt.changeTo(9, 2);

    // shadow.push(Subspace(p - 1, 1 - 1, sut, svt, 3, 3, 0.4, 0.5, 100));
    printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
    //  shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux
    while (count < 1 || !shadow.empty()) {
        count++;
        Subspace min_space = shadow.top();
        // if (min_space.mu == 1.0)
        //     break;
        shadow.pop();

        uint32_t llSize = min_space.SUSize;
        uint32_t rrSize = min_space.SVSize;
        std::vector<uint32_t> minSU(g->maxDv);
        std::vector<uint32_t> minSV(g->maxDu);
        std::unordered_map<uint32_t, uint32_t> minSUMap;
        std::unordered_map<uint32_t, uint32_t> minSVMap;
        // NOTE  remove later
        for (int i = 0; i < llSize; i++) {
            minSU[i] = min_space.SU.getByIndex(i);
        }
        for (int i = 0; i < rrSize; i++) {
            minSV[i] = min_space.SV.getByIndex(i);
        }
        ////////////////////

        for (int i = 0; i < llSize; i++) {
            minSUMap[minSU[i]] = i;
        }
        for (int i = 0; i < rrSize; i++) {
            minSVMap[minSV[i]] = i;
        }

        LinearSet candLL, candRR;
        candLL.resize(llSize + 1);
        candRR.resize(rrSize + 1);
        candLL = min_space.SU;
        candRR = min_space.SV;

        int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
        vector<int> pUU(g->maxDv * g->maxDv);
        vector<int> pVV(g->maxDu * g->maxDu);
        vector<int> e11(maxE * g->maxDv);
        vector<int> e22(maxE * g->maxDu);
        vector<uint32_t> candLtemp2(g->maxDv * g->maxDv);
        vector<uint32_t> outPosition2(g->maxDu * g->maxDu);
        // printf("====================================\n");
        vector<int> twoHopCount2(g->maxDu * g->maxDu);
        // printf("min_space.ZCount: %f\n", min_space.ZCount);
        // printf("min_space.mu: %f\n", min_space.mu);

        pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
        // printf("totalZinShadow before substracting: %d\n", totalZinShadow);
        totalZinShadow -= min_space.ZCount;
        // printf("min_space.ZCount: %f\n", min_space.ZCount);
        // printf("totalZinShadow after substracting: %f\n", totalZinShadow);
        // printf("min_space.SUSize: %d\n", min_space.SUSize);
        // printf("min_space.SVSize: %d\n", min_space.SVSize);
        // printf("min_space.p_prime: %d\n", min_space.p_prime);
        // printf("min_space.q_prime: %d\n", min_space.q_prime);
        // printf("min_space.mu: %f\n", min_space.mu);
        // printf("====================================\n");

        // if (min_space.SUSize < min_space.p_prime ||
        // 	min_space.SVSize < min_space.q_prime)
        // 	continue;

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
        e11.resize(edgeIndex);

        // Convert pVV to prefix sum array
        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }

        // Filling e22
        edgeIndex = 0;
        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                // todo use hashmap to get the index of v

                int index = minSVMap[v];
                e22[pVV[index]] = candLL[u];
                pVV[index]++;
                edgeIndex++;
            }
        }

        // Reset pVV to its prefix sum form
        // for (int i = 0; i < rrSize + 1; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        for (int v = rrSize + 1; v >= 1; v--) {
            pVV[v] = pVV[v - 1];
        }
        // for (int i = 0; i < rrSize + 1; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        pVV[0] = 0;

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////
        // printing subgraph
        // for (int i = 0; i < llSize + 1; i++) {
        //     printf("pUU[%d]: %d\n", i, pUU[i]);
        // }
        // for (int i = 0; i < llSize; i++) {
        //     for (int j = pUU[i]; j < pUU[i + 1]; j++) {
        //         printf("e11[%d]: %d\n", j, e11[j]);
        //     }
        // }
        // for (int i = 0; i < rrSize + 1; i++) {
        //     printf("pVV[%d]: %d\n", i, pVV[i]);
        // }
        // for (int i = 0; i < rrSize + 1; i++) {
        //     for (int j = pVV[i]; j < pVV[i + 1]; j++) {
        //         printf("e22[%d]: %d\n", j, e22[j]);
        //     }
        // }

        candLL.clear();
        candRR.clear();
        // todo remove 10 and add 1
        candLL.resize(g->n1 + 10);
        candRR.resize(g->n2 + 10);

        for (uint32_t uIndex = 0; uIndex < llSize; uIndex++) {
            uint32_t u = minSU[uIndex];
            printf("u: %d\n", u);
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
                    if (candLL.idx(w) < countCandLtemp2) {
                        candLL.changeTo(w, candLSize2++);
                    }
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
                // for (int i = 0; i < candLSize2; i++) {
                //     printf("candL[%d]: %d\n", i, candL[i]);
                // }

                // for (int i = 0; i < candRSize2; i++) {
                //     printf("candR[%d]: %d\n", i, candR[i]);
                // }

                int maxPossibleLength = buildDP(candLSize2, candRSize2, min_space.p_prime - 1, min_space.q_prime - 1);

                if (maxPossibleLength < minPQ2 - 1)
                    continue;

                double zInMinSubSpace = 0;

                for (uint32_t i = 0; i < pU[candLSize2]; i++) {
                    zInMinSubSpace += dpU[minPQ2 - 1][i];
                }

                totalZinShadow += zInMinSubSpace;
                if (zInMinSubSpace < 1)
                    continue;

                // sampling
                uint32_t auxSampleSize = 100;
                uint32_t auxSampleSizeTemp = auxSampleSize;
                uint32_t pqBicluqeContainZ = 0;

                while (auxSampleSizeTemp--) {
                    // selecting the first edge
                    selectedL.clear();
                    selectedR.clear();
                    double point = distribution(eng);
                    double temp = 0.0;
                    uint32_t preU = 0, preV = 0, preE = 0;
                    bool selectedLEmpty = true;

                    for (int u = 0; u < candLSize2; u++) {
                        for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                            temp += dpU[minPQ2 - 1][mapUtoV[i]];
                            if (temp + 1e-8 >= point * zInMinSubSpace) {
                                selectedR.push_back(e1[i]);
                                selectedL.push_back(u);
                                preU = u;
                                preV = e1[i];
                                preE = i;
                                selectedLEmpty = false;
                                break;
                            }
                        }
                        if (!selectedLEmpty) {
                            break;
                        }
                    }

                    assert(selectedR.size() > 0);
                    assert(selectedL.size() > 0);
                    // selecting the rest of the edge
                    for (int i = 1; i < minPQ2 - 1; i++) {
                        // TODO update to use binary search to find the vertex that start the
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

                    // end pasted
                    if (selectedR.size() < minPQ2 - 1 || selectedL.size() < minPQ2 - 1)
                        printf("selectedR.size(): %d, selectedL.size(): %d\n", selectedR.size(), selectedL.size());

                    // assert(selectedR.size() == minPQ2 - 1);
                    // assert(selectedL.size() == minPQ2 - 1);

                    // check whether the sampled z is a biclique
                    bool connected = true;
                    for (int i = 0; i < minPQ2 - 1; i++) {
                        for (int j = 0; j < minPQ2 - 1; j++) {
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
                    uint32_t noOfVsConnected = 0;

                    for (int i = preV + 1; i < candRSize2; i++) {
                        uint32_t v = candRR[i];
                        bool check = true;
                        for (int i = 0; i < minPQ2 - 1; i++) {
                            if (!g->connectUV(candLL[selectedL[i]], v)) {
                                check = false;
                                break;
                            }
                        }
                        if (check)
                            noOfVsConnected++;
                    }

                    if (noOfVsConnected >= min_space.q_prime - min_space.p_prime) {
                        pqBicluqeContainZ +=
                            C[noOfVsConnected][min_space.q_prime - min_space.p_prime];
                    }
                }
                shadow.push(Subspace(
                    min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR,
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

    printf("---------------------\n");
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
            candL = s.SU;
            candR = s.SV;
            int maxLength = buildDP(s.SUSize, s.SVSize, s.p_prime, s.q_prime);
            uint32_t llSize = s.SUSize;
            uint32_t rrSize = s.SVSize;
            LinearSet candLL, candRR;
            candLL.resize(llSize + 1);
            candRR.resize(rrSize + 1);
            candLL = s.SU;
            candRR = s.SV;
            int minPQ2 = std::min(s.p_prime, s.q_prime);
            uint zInSubSpace = s.ZCount;

            edges.clear();
            probabilitiesForFirstE.clear();
            probForFirstE.clear();
            aliasForFirstE.clear();

            for (int u = 0; u < s.SUSize; u++) {
                for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                    edge.first = u;
                    edge.second = i;
                    edges.push_back(edge);
                    probabilitiesForFirstE.push_back(dpU[minPQ2][mapUtoV[i]] / zInSubSpace);
                }
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
                selectedL.push_back(edge.first);
                selectedR.push_back(e1[edge.second]);
                preU = edge.first;
                preV = e1[edge.second];
                preE = edge.second;

                // for (int u = 0; u < s.SUSize; u++) {
                // 	for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                // 		temp += dpU[minPQ2][mapUtoV[i]];
                // 		if (temp + 1e-8 >= point * zInSubSpace) {
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

                assert(selectedR.size() == minPQ2);
                assert(selectedL.size() == minPQ2);
                // check whether the sampled z is a biclique
                bool connected = true;
                for (int i = 0; i < minPQ2; i++) {
                    for (int j = 0; j < minPQ2; j++) {
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

                uint32_t noOfVsConnected = 0;

                if (connected && p == q) {
                    noOfVsConnected = 1;
                }

                if (connected && p != q) {
                    for (int i = preV + 1; i < s.SVSize; i++) {
                        uint32_t v = candR[i];
                        bool check = true;
                        for (int i = 0; i < minPQ2; i++) {
                            if (!g->connectUV(candL[selectedL[i]], v)) {
                                check = false;
                                break;
                            }
                        }
                        if (check)
                            noOfVsConnected++;
                    }
                }
                if (noOfVsConnected >= q - p) {
                    pqBicluqeContainZ += C[noOfVsConnected][q - p];
                    successSamples++;
                }

                if (successSamples >= int(gamma)) {
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
    printf("density of biclique: %f\n", pqBicluqeContainZ / (double)t);
    printf("biclique count: %f\n", pqBicluqeContainZ / (double)t * totalZinShadow);
    printf("time: %f\n", (double)(std::chrono::high_resolution_clock::now() - start).count() / 1e9);
}

void accuracy::shadowBuilderZStar(int p, int q, double e) {
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

    // TODO use fibonacci heap
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
        // for(int i=0; i<k; i++){
        // 	printf("--------------------------------\n");
        // 	for(int j=0; j<pV[rSize]; j++){
        // 		printf("dpV[%d][%d]: %f\n", i, j, dpV[i][j]);
        // 	}
        // }
        // printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        // printf("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
        return k;
    };

    for (uint32_t u = 0; u < g->n1; u++) {
        if (g->deg1(u) < 2) {
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
            double z2 = 0;

            if (minPQ - 1 == 1) {
                for (uint32_t i = 0; i < pV[candRSize]; i++) {
                    zInSubSpace += dpV[minPQ - 1][i];
                }

            }

            else {
                for (uint32_t i = 0; i < pU[candLSize]; i++) {
                    zInSubSpace += dpU[minPQ - 1][i];
                }
            }

            totalZinShadow += zInSubSpace;

            // TODO check zInSubspace size larger than 0 and take action accordingly
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
                // 	for (int v = 0; v < candRSize; v++) {
                // 		for (int i = pV[v]; i < pV[v + 1]; i++) {
                // 			temp += dpV[minPQ - 1][mapVtoU[i]];
                // 			if (temp + 1e-8 >= point * zInSubSpace) {
                // 				selectedR.push_back(v);
                // 				selectedL.push_back(e2[i]);
                // 				preU = e2[i];
                // 				preV = v;
                // 				preE = i;
                // 				selectedLEmpty = false;
                // 				break;
                // 			}
                // 		}
                // 		if(!selectedLEmpty){
                // 			break;
                // 		}
                // 	}
                // }
                // else {
                // 	double temp = 0.0;
                // 	for (int u = 0; u < candLSize; u++) {
                // 		for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                // 			temp += dpU[minPQ - 1][mapUtoV[i]];
                // 			if (temp + 1e-8 >= point * zInSubSpace) {
                // 				selectedR.push_back(e1[i]);
                // 				selectedL.push_back(u);
                // 				preU = u;
                // 				preV = e1[i];
                // 				preE = i;
                // 				selectedLEmpty = false;
                // 				break;
                // 			}
                // 		}
                // 		if (!selectedLEmpty) {
                // 			break;
                // 		}
                // 	}
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
                for (int i = 1; i < minPQ - 1; i++) {
                    // TODO update to use binary search to find the vertex that start the
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
                        temp = 0.0;
                        point = distribution(eng);
                        for (uint32_t j = pU[preU]; j < pU[preU + 1]; j++) {
                            if (e1[j] > preV) {
                                temp += dpU[minPQ - i - 1][mapUtoV[j]];
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

                assert(index >= pU[preU] && index < pU[preU + 1]);

                std::vector<uint32_t> selectedRStar = sampeleRest(vCount, preU, pU, index);

                for (int i = 0; i < selectedRStar.size(); i++) {
                    selectedR.push_back(e1[selectedRStar[i]]);
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
            if (zInSubSpace > 0.5) {
                shadow.push(Subspace(p - 1, q - 1, candL, candR, candLSize, candRSize,
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

    // getting k paths for each edge in G is done

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
    while (shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux) {
        count++;
        Subspace min_space = shadow.top();
        // if (min_space.mu == 1.0)
        // 	break;
        shadow.pop();

        uint32_t llSize = min_space.SUSize;
        uint32_t rrSize = min_space.SVSize;
        LinearSet candLL, candRR;
        candLL.resize(llSize + 1);
        candRR.resize(rrSize + 1);
        candLL = min_space.SU;
        candRR = min_space.SV;
        int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
        vector<int> pUU(g->maxDv * g->maxDv);
        vector<int> pVV(g->maxDu * g->maxDu);
        vector<int> e11(maxE * g->maxDv);
        vector<int> e22(maxE * g->maxDu);
        vector<uint32_t> candLtemp2(g->maxDv * g->maxDv);
        vector<uint32_t> outPosition2(g->maxDu * g->maxDu);
        std::vector<uint32_t> minSU(g->maxDv);
        std::vector<uint32_t> minSV(g->maxDu);
        std::unordered_map<uint32_t, uint32_t> minSUMap;
        std::unordered_map<uint32_t, uint32_t> minSVMap;
        vector<int> twoHopCount2(g->maxDu * g->maxDu);

        pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
        totalZinShadow -= min_space.ZCount;

        // NOTE shadow has vectors instead of lineaset this can be removed
        for (int i = 0; i < llSize; i++) {
            minSU[i] = min_space.SU.getByIndex(i);
        }
        for (int i = 0; i < rrSize; i++) {
            minSV[i] = min_space.SV.getByIndex(i);
        }

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
        e11.resize(edgeIndex);

        // Convert pVV to prefix sum array
        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }

        // Filling e22
        edgeIndex = 0;
        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                // todo use hashmap to get the index of v

                int index = minSVMap[v];
                e22[pVV[index]] = candLL[u];
                pVV[index]++;
                edgeIndex++;
            }
        }

        // Reset pVV to its prefix sum form
        for (int v = rrSize + 1; v >= 1; v--) {
            pVV[v] = pVV[v - 1];
        }
        pVV[0] = 0;

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////

        candLL.clear();
        candRR.clear();
        // todo remove 10 and add 1
        candLL.resize(g->n1 + 1);
        candRR.resize(g->n2 + 1);

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
                    if (candLL.idx(w) < countCandLtemp2) {
                        candLL.changeTo(w, candLSize2++);
                    }
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

                int maxPossibleLength = buildDP(candLSize2, candRSize2, min_space.p_prime - 1, min_space.q_prime - 1);

                // printf("maxPossibleLength: %d\n", maxPossibleLength);

                if (maxPossibleLength <= minPQ2 - 1) continue;

                double zInMinSubSpace = 0;

                if (minPQ2 - 1 == 1) {
                    for (uint32_t i = 0; i < pV[candRSize2]; i++) {
                        zInMinSubSpace += dpV[minPQ2 - 1][i];
                    }

                }

                else {
                    for (uint32_t i = 0; i < pU[candLSize2]; i++) {
                        zInMinSubSpace += dpU[minPQ2 - 1][i];
                    }
                }

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
                        // TODO update to use binary search to find the vertex that start the
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

                shadow.push(Subspace(
                    min_space.p_prime - 1, min_space.q_prime - 1, candLL, candRR,
                    candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize,
                    zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize),
                    zInMinSubSpace));

                pqBicluqeEstAux +=
                    zInMinSubSpace * pqBicluqeContainZ / (double)auxSampleSize;
            }
        }
        // printf("pqBicluqeEstAux: %f\n", pqBicluqeEstAux);
        // printf("totalZinShadow: %d\n", totalZinShadow);
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
            candL = s.SU;
            candR = s.SV;

            int maxLength = buildDP(s.SUSize, s.SVSize, s.p_prime, s.q_prime);
            uint32_t llSize = s.SUSize;
            uint32_t rrSize = s.SVSize;
            LinearSet candLL, candRR;
            candLL.resize(llSize + 1);
            candRR.resize(rrSize + 1);
            candLL = s.SU;
            candRR = s.SV;
            int minPQ2 = std::min(s.p_prime, s.q_prime);
            double zInSubSpace = s.ZCount;

            edges.clear();
            probabilitiesForFirstE.clear();
            probForFirstE.clear();
            aliasForFirstE.clear();
            // fixme
            double zTemp = 0;

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
    std::uniform_int_distribution<> dis(outIndex + 1, pU[preU + 1] - 1);
    uint32_t selectedIndicesSize = 0;

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
