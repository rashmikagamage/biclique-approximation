#include "BCTV2.h"

#include <signal.h>

#include <algorithm>
#include <chrono>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../tools/hopstotchHash.hpp"
int maxPcount = 0;
int maxQcount = 0;
int maxDepth = 0;
std::vector<std::vector<uint32_t>> SUStorage;
std::vector<std::vector<uint32_t>> SVStorage;
std::vector<uint32_t> SVnew;
std::vector<uint32_t> SUnew;
void BCTV2::buildTree() {
    int maxDepth = 2 * std::max(g->maxDu, g->maxDv);
    printf("maxDu: %d\n", g->maxDu);
    printf("maxDv: %d\n", g->maxDv);
    printf("n1: %d\n", g->n1);
    printf("n2: %d\n", g->n2);
    fflush(stdout);
    SUnew.resize(g->maxDv + 1);
    SVnew.resize(g->maxDu + 1);
    SUStorage.resize(maxDepth + 1);
    SVStorage.resize(maxDepth + 1);
    int maxVec = 2 * std::max(g->maxDu, g->maxDv);
    count.resize(maxVec + 1, std::vector<double>(maxVec + 1, 0));
    maxDepth = 2 * std::max(g->maxDu, g->maxDv);
    NodeV root;
    auto start = std::chrono::high_resolution_clock::now();

    std::vector<uint32_t> SU;
    std::vector<uint32_t> SV;
    // std::vector<uint32_t> twoHopCountU(g->n1 + 1, 0);
    // std::vector<uint32_t> twoHopCountV(g->n2 + 1, 0);
    std::vector<uint32_t> twoHopCount(g->n1 + 1, 0);
    // std::vector<uint32_t> twoHopTempU;
    // std::vector<uint32_t> twoHopTempV;
    std::vector<uint32_t> twoHopTemp;
    int maxUProcessed = -1;
    int maxVProcessed = -1;

    for (uint32_t u = 0; u < g->n1; u++) {
        SU.clear();
        SV.clear();

        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            SV.push_back(v);

            for (uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if (w == u) break;
                if (twoHopCount[w]++ == 0) {
                    SU.push_back(w);
                    twoHopTemp.push_back(w);
                }
            }
        }
        for (int i = 0; i < twoHopTemp.size(); i++) {
            twoHopCount[twoHopTemp[i]] = 0;
            twoHopTemp[i] = 0;
        }

        if (!SU.empty() && !SV.empty()) {
            NodeV* n = new NodeV(2, 1, u);
            SUStorage[0] = SU;
            SVStorage[0] = SV;
            n->uh = 1;
            n->depth = 0;
            processNode(n);
        }
    }
    // for (uint16_t i = 0; i < g->n1 + g->n2; i++) {
    //     SU.clear();
    //     SV.clear();
    //     int vertex = g->degeneracy[i];
    //     int side = g->side[i];
    //     if (side == 1) {
    //         maxUProcessed = vertex;
    //         for (uint32_t j = g->pU[vertex]; j < g->pU[vertex + 1]; j++) {
    //             uint32_t v = g->e1[j];
    //             if (v <= maxVProcessed) continue;
    //             SV.push_back(v);
    //             for (uint32_t k = g->pV[v + 1]; k > g->pV[v]; k--) {
    //                 uint32_t w = g->e2[k - 1];
    //                 if (w == vertex) break;
    //                 if (twoHopCountU[w]++ == 0) {
    //                     SU.push_back(w);
    //                     twoHopTempU.push_back(w);
    //                 }
    //             }
    //         }
    //         for (int i = 0; i < twoHopTempU.size(); i++) {
    //             twoHopCountU[twoHopTempU[i]] = 0;
    //             twoHopTempU[i] = 0;
    //         }
    //         if (SU.empty() && SV.empty()) continue;
    //         NodeV* n = new NodeV(SU, SV, 2, 1, vertex);
    //         n->uh = 1;
    //         processNode(n);
    //     } else {
    //         maxVProcessed = vertex;
    //         for (uint32_t j = g->pV[vertex]; j < g->pV[vertex + 1]; j++) {
    //             uint32_t u = g->e2[j];
    //             if (u <= maxUProcessed) continue;
    //             SU.push_back(u);
    //             for (uint32_t k = g->pU[u + 1]; k > g->pU[u]; k--) {
    //                 uint32_t w = g->e1[k - 1];
    //                 if (w == vertex) break;
    //                 if (twoHopCountV[w]++ == 0) {
    //                     SV.push_back(w);
    //                     twoHopTempV.push_back(w);
    //                 }
    //             }
    //         }
    //         for (int i = 0; i < twoHopTempV.size(); i++) {
    //             twoHopCountV[twoHopTempV[i]] = 0;
    //             twoHopTempV[i] = 0;
    //         }
    //         if (SU.empty() && SV.empty()) continue;
    //         NodeV* n = new NodeV(SU, SV, 2, 2, vertex);
    //         n->vh = 1;
    //         processNode(n);
    //     }
    // }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Tree completed in: %ld (ms)\n", duration.count() / 1000);

    printf("count size: %d\n", count.size());
    // std::vector<std::vector<double>> result = countBicliques(&root);
    printf("maxPcount: %d\n", maxPcount);
    printf("maxQcount: %d\n", maxQcount);
    for (int i = 2; i <= maxPcount; i++) {
        for (int j = 2; j <= maxQcount; j++) {
            if (count[i][j] > 0.5) printf("(%d - %d) = %f\n", i, j, count[i][j]);
        }
    }
}

void BCTV2::processNode(NodeV* node) {
    int depth = node->depth;

    if (SUStorage[depth].empty() && SVStorage[depth].empty()) {
        int uh = node->uh;
        int vh = node->vh;
        int up = node->up;
        int vp = node->vp;

        for (int i = 0; i <= up; ++i) {
            for (int j = 0; j <= vp; ++j) {
                maxPcount = std::max(maxPcount, i + uh);
                maxQcount = std::max(maxQcount, j + vh);
                int maxSize = std::max(g->maxDu, g->maxDv) * 2;
                // assert(i + uh < maxSize);
                // assert(j + vh < maxSize);
                // assert(up < maxSize);
                // assert(vp < maxSize);
                count[i + uh][j + vh] += C[up][i] * C[vp][j];
            }
        }
        return;
    }
    // printf("node->SU size: %d\n", node->SU.size());
    if (!SUStorage[depth].empty() && !SVStorage[depth].empty()) {
        std::pair<uint32_t, bool> pivotPair = selectPivot(SUStorage[depth], SVStorage[depth]);
        uint32_t pivot = pivotPair.first;

        if (pivotPair.second) {
            std::vector<uint32_t> SVPrime;
            // std::unordered_set<uint32_t> SUSet;
            // std::unordered_set<uint32_t> SVSet;
            hopstotchHash SUSet;
            hopstotchHash SVSet;
            for (uint32_t v : SVStorage[depth]) {
                SVSet.insert(v);
                if (g->connectUV(pivot, v)) {
                    SVnew.push_back(v);

                } else {
                    SVPrime.push_back(v);
                }
            }

            for (uint32_t u : SUStorage[depth]) {
                SUSet.insert(u);

                if (u != pivot) {
                    SUnew.push_back(u);
                }
            }

            NodeV* newNode = new NodeV(1, 1, pivot);
            SUStorage[depth + 1] = std::move(SUnew);
            SVStorage[depth + 1] = std::move(SVnew);
            newNode->depth = node->depth + 1;
            newNode->up += node->up + 1;
            newNode->vh = node->vh;
            newNode->uh = node->uh;
            newNode->vp = node->vp;

            processNode(newNode);

            hopstotchHash visited;

            for (uint32_t v_i : SVPrime) {
                visited.insert(v_i);
                std::vector<uint32_t> SVchild;
                std::vector<uint32_t> SUchild;
                std::vector<uint32_t> twoHopsCount(g->n2 + 1);
                std::unordered_set<uint32_t> twoHops;

                // /hopstotchHash twoHopsCount
                for (uint32_t j = g->pV[v_i]; j < g->pV[v_i + 1]; j++) {
                    uint32_t u = g->e2[j];

                    if (SUSet.contain(u) == 1) {
                        SUchild.push_back(u);
                    }

                    for (uint32_t k = g->pU[u]; k < g->pU[u + 1]; k++) {
                        uint32_t v = g->e1[k];

                        if (SVSet.contain(v) == 1 && visited.contain(v) == 0 && twoHopsCount[v] == 0) {
                            SVchild.push_back(v);
                        }
                        twoHopsCount[v]++;
                    }
                }

                // if (g->pV[v_i + 1] - g->pV[v_i] < SU.size()) {
                // } else {
                //     for (uint32_t u : SU) {
                //         if (g->connectUV(u, v_i)) {
                //             SUchild.push_back(u);
                //         }
                //     }
                //     for (int j = g->pV[v_i]; j < g->pV[v_i + 1]; j++) {
                //         uint32_t u = g->e2[j];
                //         for (int k = g->pU[u]; k < g->pU[u + 1]; k++) {
                //             uint32_t v = g->e1[k];
                //             if (SVSet.contain(v) == 1 && visited.contain(v) == 0 && twoHopsCount[v] == 0) {
                //                 SVchild.push_back(v);
                //             }
                //             twoHopsCount[v]++;
                //         }
                //     }
                // }

                NodeV* childNode = new NodeV(2, 2, v_i);
                SUStorage[depth + 1] = SUchild;
                SVStorage[depth + 1] = SVchild;
                childNode->depth = node->depth + 1;
                childNode->vh += node->vh + 1;
                childNode->uh = node->uh;
                childNode->vp = node->vp;
                childNode->up = node->up;
                // node->children.push_back(childNode);
                processNode(childNode);
            }
        } else {
            std::vector<uint32_t> SUPrime;
            hopstotchHash SUSet;
            hopstotchHash SVSet;
            for (uint32_t u : SUStorage[depth]) {
                SUSet.insert(u);
                if (g->connectUV(u, pivot)) {
                    SUnew.push_back(u);
                } else {
                    SUPrime.push_back(u);
                }
            }

            for (uint32_t v : SVStorage[depth]) {
                SVSet.insert(v);
                if (v != pivot) {
                    SVnew.push_back(v);
                }
            }

            NodeV* newNode = new NodeV(1, 2, pivot);

            newNode->vp += node->vp + 1;
            newNode->vh = node->vh;
            newNode->uh = node->uh;
            newNode->up = node->up;
            SUStorage[depth + 1] = std::move(SUnew);
            SVStorage[depth + 1] = std::move(SVnew);
            newNode->depth = node->depth + 1;
            processNode(newNode);

            hopstotchHash visited;

            for (uint32_t u_i : SUPrime) {
                visited.insert(u_i);
                std::vector<uint32_t> SUchild;
                std::vector<uint32_t> SVchild;
                std::vector<uint32_t> twoHopsCount(g->n1 + 1);
                // std::unordered_set<uint16_t> twoHops;

                for (uint32_t j = g->pU[u_i]; j < g->pU[u_i + 1]; j++) {
                    uint32_t v = g->e1[j];
                    if (SVSet.contain(v) == 1) {
                        SVchild.push_back(v);
                    }
                    for (uint32_t k = g->pV[v]; k < g->pV[v + 1]; k++) {
                        uint32_t u = g->e2[k];
                        if (SUSet.contain(u) == 1 && visited.contain(u) == 0 && twoHopsCount[u] == 0) {
                            SUchild.push_back(u);
                        }
                        twoHopsCount[u]++;
                        // twoHops.insert(u);
                    }
                }

                NodeV* childNode = new NodeV(2, 1, u_i);
                // node->children.push_back(childNode);
                childNode->uh += node->uh + 1;
                childNode->vp = node->vp;
                childNode->vh = node->vh;
                childNode->up = node->up;
                SUStorage[depth + 1] = SUchild;
                SVStorage[depth + 1] = SVchild;
                childNode->depth = node->depth + 1;
                processNode(childNode);
            }
            // for (uint32_t v : SV) {
            //     SVBitMask[v] = 0;
            // }
            // for (uint32_t u : SU) {
            //     SUBitMask[u] = 0;
            // }
        }
    }

    if (SUStorage[depth].empty() || SVStorage[depth].empty()) {
        std::vector<uint32_t> SVnew;
        std::vector<uint32_t> SUnew;
        int side = 1;
        int vertex = 0;
        int pCount = 0;
        if (!SUStorage[depth].empty()) {
            pCount = SUStorage[depth].size();

        } else {
            pCount = SVStorage[depth].size();
            side = 2;
        }
        NodeV* newNode = new NodeV(side, vertex, pCount);
        if (side == 1) {
            newNode->up += node->up + pCount;
            newNode->uh = node->uh;
            newNode->vp = node->vp;
            newNode->vh = node->vh;
        } else {
            newNode->vp += node->vp + pCount;
            newNode->vh = node->vh;
            newNode->uh = node->uh;
            newNode->up = node->up;
        }
        SUStorage[depth + 1] = SUnew;
        SVStorage[depth + 1] = SVnew;
        newNode->depth = node->depth + 1;
        // node->children.push_back(newNode);
        processNode(newNode);
    }
}

std::pair<uint32_t, bool> BCTV2::selectPivot(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV) {
    auto Start = std::chrono::high_resolution_clock::now();
    uint32_t bestVertex = 0;
    bool isFromSU = true;
    uint32_t maxSum = 0;

    for (uint32_t v : SV) {
        int pU = 0;
        for (uint32_t u : SU) {
            if (g->connectUV(u, v)) {
                pU++;
            }
        }
        uint32_t sum = pU + SV.size() - 1;
        if (sum >= maxSum) {
            maxSum = sum;
            bestVertex = v;
            isFromSU = false;
        }
    }

    for (uint32_t u : SU) {
        int pV = 0;
        for (uint32_t v : SV) {
            if (g->connectUV(u, v)) {
                pV++;
            }
        }

        uint32_t sum = SU.size() - 1 + pV;
        if (sum >= maxSum) {
            maxSum = sum;
            bestVertex = u;
            isFromSU = true;
        }
    }

    auto End = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(End - Start);

    return {bestVertex, isFromSU};
}

void destroyVector(std::vector<uint32_t>& vec) {
    std::vector<uint32_t>().swap(vec);
}