#include "accuracy.h"
#include "../biGraph/biGraph.hpp"
#include <cassert>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <queue>
#include <random>
#include <unordered_map>
#include <unordered_set>
#include <vector>
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
        : p_prime(p_prime), q_prime(q_prime), SU(SU), SV(SV), SUSize(SUSize),
        SVSize(SVSize), mu(mu), countPQMean(countPQMean), ZCount(ZCount) {}
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

    std::vector<uint32_t> pV(g->maxDu); // Adjusted size to accommodate pV[v + 1]
    std::vector<uint32_t> pU(g->maxDv);
    std::vector<uint32_t> e1(maxE); // Assuming maximum possible size
    std::vector<uint32_t> e2(
        maxE); // Assuming maximum possible sizmaxPossibleLengthe
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
            }
            else {
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

        if (pU[lSize] == 0 || pV[rSize] == 0)

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
        // if no vertices found with two hop neighbors
        if (countCandLtemp == 0) {
            continue;
        }

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            int candLSize = 0;
            for (uint32_t j = outPosition[i - g->pU[u]] + 1; j < g->pV[v + 1]; j++) {
                uint32_t w = g->e2[j];
                // check whether the w is inside the candLTemo
                if (candL.idx(w) < countCandLtemp) {
                    candL.changeTo(w, candLSize++);
                }
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
            if (maxPossibleLength < minPQ - 1)
                continue;

            double zInSubSpace = 0;
            for (uint32_t i = 0; i < pU[candLSize]; i++) {
                zInSubSpace += dpU[minPQ - 1][i];
            }
            totalZinShadow += zInSubSpace;
            // TODO check zInSubspace size larger than 0 and take action accordingly
            if (zInSubSpace == 0)
                continue;
            double pqBicluqeContainZ = 0;
            double auxSampleSize = 1000;
            int auxSampleSize_temp = auxSampleSize;
            double t_sample = 0;

            auto auxSamplingStart = chrono::high_resolution_clock::now();

            while (auxSampleSize_temp--) {

                // selecting the first edge
                std::random_device rd;

                std::mt19937 eng(rd()); // Engine (Mersenne Twister)
                std::uniform_real_distribution<> distribution(
                    0.0, 1.0); // Distribution for doubles between 0.0 and 1.0
                // (inclusive of 0.0 and exclusive of 1.0)
                vector<uint32_t> selectedL(g->maxDv + 1), selectedR(g->maxDu + 1);
                selectedL.clear();
                selectedR.clear();
                double point = distribution(eng);
                double temp = 0.0;
                uint32_t preU = 0, preV = 0, preE = 0;
                bool selectedLEmpty = true;

                for (int u = 0; u < candLSize; u++) {
                    for (uint32_t i = pU[u]; i < pU[u + 1]; i++) {
                        temp += dpU[minPQ - 1][mapUtoV[i]];
                        if (temp + 1e-8 >= point * zInSubSpace) {
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
    std::mt19937 eng(rd()); // Engine (Mersenne Twister)
    std::uniform_real_distribution<> distribution(
        0.0, 1.0); // Distribution for doubles between 0.0 and 1.0 (inclusive of
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

    //shadowDuration < gamma * (tTotal / auxTotalSamples) * totalZinShadow / pqBicluqeEstAux
    int x = 2;
    while (x < 3) {
        x++;
        printf("shadowDuration: %f\n", shadowDuration);
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

        std::fill(pVV.begin(), pVV.begin() + rrSize + 1, 0);

        for (int j = 0; j < llSize; j++) {
            uint32_t x = candLL[j];
            pUU[j + 1] = pUU[j];

            for (int k = 0; k < rrSize; k++) {
                uint32_t y = candRR[k];
                if (g->connectUV(x, y)) {

                    e11[pUU[j + 1]++] = k;
                    pVV[k + 1]++;
                }
            }
        }
        if (pUU[llSize] == 0)
            continue;

        for (int v = 0; v < rrSize; v++) {
            pVV[v + 1] += pVV[v];
        }

        for (int u = 0; u < llSize; u++) {
            for (int i = pUU[u]; i < pUU[u + 1]; i++) {
                int v = e11[i];
                e22[pVV[v]] = u;
                pVV[v]++;
            }
        }
        for (int v = rrSize; v >= 1; v--) {
            pVV[v] = pVV[v - 1];
        }
        pVV[0] = 0;

        /////////////////////////////
        /* ***subgraph is done*** */
        ///////////////////////////

        for (uint32_t u = 0; u < llSize; u++) {

            uint32_t candRSize2 = 0;

            for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {
                uint32_t v = e11[i];
                candRR.changeTo(v, candRSize2);
                candRSize2++;
            }

            candLtemp2.clear();

            for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {

                uint32_t v = e11[i];
                for (uint32_t j = pVV[v + 1] - 1; j >= pVV[v]; j--) {
                    uint32_t w = e22[j];
                    if (w == u) {
                        // keepig the index of u so it can later get vertices larger than u
                        outPosition2[i - pUU[u]] = j;
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
                    candLL.changeTo(w, countCandLtemp2);
                    countCandLtemp2++;
                }
                twoHopCount2[w] = 0;
            }

            // if no vertices found with two hop neighbors

            if (countCandLtemp2 == 0) {
                continue;
            }

            for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {
                uint32_t v = e11[i];
                int candLSize2 = 0;
                uint32_t j = outPosition2[i - pUU[u]] + 1;
                for (; j < pVV[v + 1]; j++) {
                    uint32_t w = e22[j];
                    // check whether the w is inside the candL
                    if (candLL.idx(w) < countCandLtemp2) {
                        candLL.changeTo(w, candLSize2++);
                    }
                }
                candRR.changeTo(v, --candRSize2);

                if (candLSize2 < min_space.p_prime - 1)
                    continue;
                if (candRSize2 < min_space.q_prime - 1)
                    break;

                for (int i = 1; i < candRSize2; i++) {
                    candRR.changeToByPos(i, i - 1);
                }
                candL.clear();
                candR.clear();
                candL = candLL;
                candR = candRR;
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


                    //pasted
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

                    //end pasted

                    assert(selectedR.size() == minPQ2 - 1);
                    assert(selectedL.size() == minPQ2 - 1);

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
            }
            else {
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
                if (selectedL.size() != minPQ2 || selectedR.size() != minPQ2) {
                    printf("minPQ2: %d\n", minPQ2);
                    printf("selectedR: %d\n", selectedR.size());
                    printf("selectedL: %d\n", selectedL.size());
                    exit(0);
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
        }
        else {
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
        }
        else {
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
    }
    else {
        return Alias[i];
    }
}

// uint32_t llSize = min_space.SUSize;
// uint32_t rrSize = min_space.SVSize;
// LinearSet candLL, candRR;
// candLL.resize(llSize + 1);
// candRR.resize(rrSize + 1);
// candLL = min_space.SU;
// candRR = min_space.SV;
// int minPQ2 = std::min(min_space.p_prime, min_space.q_prime);
// vector<int> pUU(g->maxDv* g->maxDv);
// vector<int> pVV(g->maxDu* g->maxDu);
// vector<int> e11(maxE* g->maxDv);
// vector<int> e22(maxE* g->maxDu);
// vector<uint32_t> candLtemp2(g->maxDv* g->maxDv);
// vector<uint32_t> outPosition2(g->maxDu* g->maxDu);
// //printf("====================================\n");
// vector<int> twoHopCount2(g->maxDu* g->maxDu);
// pqBicluqeEstAux -= min_space.ZCount * min_space.mu;
// //printf("totalZinShadow before substracting: %d\n", totalZinShadow);
// totalZinShadow -= min_space.ZCount;
// // printf("min_space.ZCount: %d\n", min_space.ZCount);
// // printf("totalZinShadow after substracting: %d\n", totalZinShadow);
// // printf("min_space.SUSize: %d\n", min_space.SUSize);
// // printf("min_space.SVSize: %d\n", min_space.SVSize);
// // printf("min_space.p_prime: %d\n", min_space.p_prime);
// // printf("min_space.q_prime: %d\n", min_space.q_prime);
// // printf("min_space.mu: %f\n", min_space.mu);
// // printf("====================================\n");
// if (min_space.SUSize < min_space.p_prime || min_space.SVSize <
// min_space.q_prime) 	continue;

// /////////////////////////////
// /* ***subgraph creating*** */
// ///////////////////////////

// std::fill(pVV.begin(), pVV.begin() + rrSize + 1, 0);

// for (int j = 0; j < llSize; j++) {
// 	uint32_t x = candLL[j];
// 	pUU[j + 1] = pUU[j];

// 	for (int k = 0; k < rrSize; k++) {
// 		uint32_t y = candRR[k];
// 		if (g->connectUV(x, y)) {

// 			e11[pUU[j + 1]++] = k;
// 			pVV[k + 1]++;
// 		}
// 	}

// }
// if (pUU[llSize] == 0) continue;;

// for (int v = 0; v < rrSize; v++) {
// 	pVV[v + 1] += pVV[v];
// }

// for (int u = 0; u < llSize; u++) {
// 	for (int i = pUU[u]; i < pUU[u + 1]; i++) {
// 		int v = e11[i];
// 		e22[pVV[v]] = u;
// 		pVV[v]++;

// 	}
// }
// for (int v = rrSize; v >= 1; v--) {
// 	pVV[v] = pVV[v - 1];
// }
// pVV[0] = 0;

// /////////////////////////////
// /* ***subgraph is done*** */
// ///////////////////////////

// for (uint32_t u = 0; u < llSize; u++) {

// 	uint32_t candRSize2 = 0;

// 	for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {
// 		uint32_t v = e11[i];
// 		candRR.changeTo(v, candRSize2);
// 		candRSize2++;
// 	}

// 	candLtemp2.clear();

// 	for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {

// 		uint32_t v = e11[i];
// 		for (uint32_t j = pVV[v + 1] - 1; j >= pVV[v]; j--) {
// 			uint32_t w = e22[j];
// 			if (w == u) {
// 				// keepig the index of u so it can later get
// vertices larger than u 				outPosition2[i - pUU[u]] = j; 				break;
// 			}
// 			if (twoHopCount2[w] == 0) {
// 				candLtemp2.push_back(w);
// 			}
// 			twoHopCount2[w]++;
// 		}
// 	}

// 	int countCandLtemp2 = 0;
// 	for (auto w : candLtemp2) {
// 		if (twoHopCount2[w] > 0) {
// 			candLL.changeTo(w, countCandLtemp2);
// 			countCandLtemp2++;
// 		}
// 		twoHopCount2[w] = 0;
// 	}

// 	//if no vertices found with two hop neighbors

// 	// if (countCandLtemp2 == 0) {
// 	// 	continue;
// 	// }

// 	for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {
// 		uint32_t v = e11[i];
// 		int candLSize2 = 0;
// 		uint32_t j = outPosition2[i - pUU[u]] + 1;
// 		for (;j < pVV[v + 1]; j++) {
// 			uint32_t w = e22[j];
// 			// check whether the w is inside the candL
// 			if (candLL.idx(w) < countCandLtemp2) {
// 				candLL.changeTo(w, candLSize2++);

// 			}
// 		}
// 		candRR.changeTo(v, --candRSize2);
// 		// if(candLSize2 > 0){
// 		// 	printf("candLSize2: %d\n", candLSize2);
// 		// 	printf("candRSize2: %d\n", candRSize2);
// 		// }

// 		if (candLSize2 < min_space.p_prime - 1)
// 			continue;
// 		if (candRSize2 < min_space.q_prime - 1)
// 			break;

// 		for (int i = 1; i < candRSize2; i++) {
// 			candRR.changeToByPos(i, i - 1);
// 		}
// 		candL.clear();
// 		candR.clear();
// 		candL = candLL;
// 		candR = candRR;

// 		int maxPossibleLength = buildDP(candLSize2, candRSize2,
// min_space.p_prime - 1, min_space.q_prime - 1);

// 		if (maxPossibleLength < minPQ2 - 1) continue;

// 		uint32_t zInMinSubSpace = 0;

// 		for (uint32_t i = 0; i < pU[candLSize2]; i++)
// 		{
// 			zInMinSubSpace += dpU[minPQ2 - 1][i];
// 		}

// 		totalZinShadow += zInMinSubSpace;

// 		// printf("candLSize2: %d\n", candLSize2);
// 		// printf("candRSize2: %d\n", candRSize2);
// 		// printf("min_space.p_prime: %d\n", min_space.p_prime);
// 		// printf("min_space.q_prime: %d\n", min_space.q_prime);
// 		// printf("zInMinSubSpace: %d\n", zInMinSubSpace);
// 		// printf("totalZinShadow after adding from min space: %d\n",
// totalZinShadow);

// 		if (zInMinSubSpace == 0) continue;

// 		//sampling
// 		uint32_t auxSampleSize = 10;
// 		uint32_t auxSampleSizeTemp = auxSampleSize;
// 		uint32_t pqBicluqeContainZ = 0;
// 		while (auxSampleSizeTemp--) {
// 			//selecting the first edge
// 			selectedL.clear();
// 			selectedR.clear();
// 			double point = distribution(eng);
// 			double temp = 0.0;
// 			uint32_t preU = 0, preV = 0, preE = 0;
// 			bool selectedLEmpty = true;
// 			for (int u = 0; u < candLSize2; u++) {
// 				for (uint32_t i = pUU[u]; i < pUU[u + 1]; i++) {
// 					temp += dpU[minPQ2 - 1][mapUtoV[i]];
// 					if (temp + 1e-8 >= point *
// zInMinSubSpace) { 						selectedR.push_back(e1[i]); 						selectedL.push_back(u); 						preU =
// u; 						preV = e1[i]; 						preE = i; 						selectedLEmpty = false; 						break;
// 					}
// 				}
// 				if (!selectedLEmpty) {
// 					break;
// 				}
// 			}

// 			//assert(selectedR.size() > 0);
// 			//selecting the rest of the edge
// 			for (int i = 1; i < minPQ2 - 1; i++) {

//
// 				//select u
// 				temp = 0.0;
// 				point = distribution(eng);
// 				for (uint32_t j = pVV[preV]; j < pVV[preV + 1];
// j++)
// 				{
// 					if (e22[j] > preU) {
// 						temp += dpV[minPQ2 - i -
// 1][mapVtoU[j]]; 						if (temp + 1e-8 >= point * dpU[minPQ - i][mapUtoV[preE]]) {
// 							uint32_t u = e22[j];
// 							selectedL.push_back(u);
// 							preU = u;
// 							preE = j;
// 							break;
// 						}
// 					}
// 				}
// 				//select v
// 				temp = 0.0;
// 				point = distribution(eng);
// 				for (uint32_t j = pUU[preU]; j < pUU[preU + 1];
// j++) { 					if (e11[j] > preV) { 						temp += dpU[minPQ2 - i - 1][mapUtoV[j]]; 						if (temp
// + 1e-8 >= point * dpV[minPQ - i - 1][mapVtoU[preE]]) {

// 							uint32_t v = e11[j];
// 							selectedR.push_back(v);
// 							preV = v;
// 							preE = j;
// 							break;
// 						}
// 					}
// 				}
// 			}

// 			if (selectedL.size() < minPQ2 - 1) continue;
// 			//assert(selectedR.size() == minPQ2 - 1);
// 			//assert(selectedL.size() == minPQ2 - 1);

// 			//check whether the sampled z is a biclique
// 			bool connected = true;
// 			for (int i = 0; i < minPQ2 - 1; i++) {
// 				for (int j = 0; j < minPQ2 - 1; j++) {
// 					if (i == j) continue;
// 					if (i > 0 && i == j + 1) continue;

// 					if (!g->connectUV(candL[selectedL[i]],
// candR[selectedR[j]])) { 						connected = false; 						break;
// 					}
// 				}
// 				if (!connected) break;
// 			}

// 			if (!connected) continue;
// 			uint32_t noOfVsConnected = 0;

// 			for (int i = preV + 1; i < candRSize2; i++) {
// 				uint32_t v = candRR[i];
// 				bool check = true;
// 				for (int i = 0; i < minPQ2 - 1; i++) {
// 					if (!g->connectUV(candLL[selectedL[i]],
// v)) { 						check = false; 						break;
// 					}
// 				}
// 				if (check) noOfVsConnected++;

// 			}

// 			if (noOfVsConnected >= min_space.q_prime -
// min_space.p_prime) { 				pqBicluqeContainZ +=
// C[noOfVsConnected][min_space.q_prime - min_space.p_prime];
// 			}

// 		}
// 		shadow.push(Subspace(min_space.p_prime - 1, min_space.q_prime -
// 1, candLL, candRR, candLSize2, candRSize2, pqBicluqeContainZ / auxSampleSize,
// zInMinSubSpace * (pqBicluqeContainZ / auxSampleSize), zInMinSubSpace));
// 		pqBicluqeEstAux += zInMinSubSpace * pqBicluqeContainZ /
// (double)auxSampleSize;

// 	}
// }

// void initializeAliasMethod(const std::vector<double>& weights,
// std::vector<int>& alias, std::vector<double>& prob) { 	int n = weights.size();
// 	std::vector<double> normalized_weights(weights);
// 	double sum = 0;
// 	for (double w : weights) sum += w;
// 	for (double& w : normalized_weights) w *= n / sum;

// 	std::vector<int> small, large;
// 	for (int i = 0; i < n; ++i) {
// 		if (normalized_weights[i] < 1.0)
// 			small.push_back(i);
// 		else
// 			large.push_back(i);
// 	}

// 	alias.resize(n);
// 	prob.resize(n);

// 	while (!small.empty() && !large.empty()) {
// 		int l = small.back(), g = large.back();
// 		small.pop_back();
// 		large.pop_back();

// 		prob[l] = normalized_weights[l];
// 		alias[l] = g;

// 		normalized_weights[g] = (normalized_weights[g] +
// normalized_weights[l]) - 1; 		if (normalized_weights[g] < 1.0)
// 			small.push_back(g);
// 		else
// 			large.push_back(g);
// 	}

// 	while (!large.empty()) {
// 		int g = large.back();
// 		large.pop_back();
// 		prob[g] = 1.0;
// 	}

// 	while (!small.empty()) {
// 		int l = small.back();
// 		small.pop_back();
// 		prob[l] = 1.0;
// 	}
// }

// int sample(const std::vector<int>& alias, const std::vector<double>& prob) {
// 	int n = prob.size();
// 	std::random_device rd;
// 	std::mt19937 gen(rd());
// 	std::uniform_int_distribution<> dis(0, n - 1);
// 	std::uniform_real_distribution<> dis_real(0.0, 1.0);

// 	int i = dis(gen);
// 	double r = dis_real(gen);
// 	return r < prob[i] ? i : alias[i];
// }
/*
auto computeDP = [&](int leftSize, int rightSize) -> int {
        std::fill(pV.begin(), pV.begin() + rightSize + 1, 0);
        for (int leftIndex = 0; leftIndex < leftSize; leftIndex++) {
                uint32_t leftNode = candL[leftIndex];
                pU[leftIndex + 1] = pU[leftIndex];

                if (rightSize < g->deg1(leftNode)) {
                        for (int rightIndex = 0; rightIndex < rightSize;
rightIndex++) { uint32_t rightNode = candR[rightIndex]; if
(g->connectUV(leftNode, rightNode)) { e1[pU[leftIndex + 1]++] = rightIndex;
                                        pV[rightIndex + 1]++;
                                }
                        }
                }
                else {
                        auto start = g->e1.begin() + g->pU[leftNode];
                        auto end = g->e1.begin() + g->pU[leftNode + 1];
                        uint32_t i = std::upper_bound(start, end, candR[0]) -
g->e1.begin(); for (; i < g->pU[leftNode + 1]; i++) { int rightIndex =
candR.idx(g->e1[i]); if (rightIndex < rightSize) { e1[pU[leftIndex + 1]++] =
rightIndex; pV[rightIndex + 1]++;
                                }
                        }
                }
        }

        if (pU[leftSize] == 0) return 0;

        for (int rightIndex = 0; rightIndex < rightSize; rightIndex++) {
                pV[rightIndex + 1] += pV[rightIndex];
        }

        for (int leftIndex = 0; leftIndex < leftSize; leftIndex++) {
                for (int i = pU[leftIndex]; i < pU[leftIndex + 1]; i++) {
                        int rightIndex = e1[i];
                        e2[pV[rightIndex]] = leftIndex;
                        mapUtoV[i] = pV[rightIndex];
                        mapVtoU[pV[rightIndex]] = i;
                        pV[rightIndex]++;
                }
        }

        std::fill(pV.begin(), pV.begin() + rightSize + 1, 0);

        int minLR = std::min(leftSize, rightSize);
        int k = 2;
        for (; k <= minLR && k < minPQ; k++) {
                std::fill(ddp.begin(), ddp.begin() + pU[leftSize], 0);
                updateDPValues(rightSize, k, dpV, dpU, ddp, pV, e2, mapVtoU);
                std::fill(ddp.begin(), ddp.begin() + pU[leftSize], 0);
                updateDPValues(leftSize, k, dpU, dpV, ddp, pU, e1, mapUtoV);
        }

        return k;
        };

void updateDPValues(int size, int k, std::vector<double>& sourceDP,
std::vector<double>& targetDP, std::vector<double>& ddp, std::vector<int>&
positions, std::vector<int>& edges, std::vector<int>& map) { for (int index = 0;
index < size; index++) { for (int i = positions[index] + 1; i < positions[index
+ 1]; i++) { int node = edges[i]; if (sourceDP[k - 1][map[i]] < 0.5) continue;
                        ddp[positions[index]] += sourceDP[k - 1][map[i]];
                        ddp[i] -= sourceDP[k - 1][map[i]];
                }
                for (int e = positions[index]; e < positions[index + 1]; e++) {
                        targetDP[k][e] = std::max(0.0, targetDP[k][e - 1] +
ddp[e]);
                }
        }
}

*/