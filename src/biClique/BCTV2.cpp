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

void BCTV2::buildTree() {
    int maxVec = std::max(g->maxDu, g->maxDv);
    count = std::vector<std::vector<double>>(100, std::vector<double>(100, 0));

    // SUBitMask = std::vector<uint32_t>(g->n1 + 2, 0);
    // SVBitMask = std::vector<uint32_t>(g->n2 + 2, 0);

    NodeV root;
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<uint32_t> SU(g->n1 + 1);
    std::vector<uint32_t> SV(g->n2 + 1);
    // std::stack<NodeV*> Q;
    std::vector<uint32_t> twoHopCount(g->n1 + 1, 0);

    for (uint32_t u = 0; u < g->n1; u++) {
        SU.clear();
        SV.clear();

        std::vector<uint32_t> twoHopCount(g->n1 + 1);
        for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
            uint32_t v = g->e1[i];
            SV.push_back(v);

            for (uint32_t j = g->pV[v + 1] - 1; j >= g->pV[v]; j--) {
                uint32_t w = g->e2[j];
                if (w == u) {
                    break;
                }
                if (twoHopCount[w] == 0) {
                    SU.push_back(w);
                }
                twoHopCount[w]++;
            }
        }

        NodeV* n = new NodeV(SU, SV, 2, 1, u);
        SU.shrink_to_fit();
        SV.shrink_to_fit();
        n->uh = 1;

        processNode(n);
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    printf("Tree completed in: %ld (ms)\n", duration.count() / 1000);

    // std::vector<std::vector<double>> result = countBicliques(&root);
    for (int i = 2; i < 5; i++) {
        for (int j = 2; j < 5; j++) {
            printf("(%d - %d) = %f\n", i, j, count[i][j]);
        }
    }
}

void BCTV2::processNode(NodeV* node) {
    std::vector<uint32_t> SU = node->SU;
    std::vector<uint32_t> SV = node->SV;

    if (SU.empty() && SV.empty()) {
        int uh = node->uh;
        int vh = node->vh;
        int up = node->up;
        int vp = node->vp;
        int upCap = std::min(100, up);
        int vpCap = std::min(100, vp);
        if (uh > 100 || vh > 100) return;
        for (int i = 0; i <= upCap; ++i) {
            for (int j = 0; j <= vpCap; ++j) {
                count[i + uh][j + vh] += C[up][i] * C[vp][j];
            }
        }
        return;
    }

    if (!SV.empty() && !SU.empty()) {
        std::pair<uint32_t, bool> pivotPair = selectPivot(SU, SV);
        uint32_t pivot = pivotPair.first;

        if (pivotPair.second) {
            std::vector<uint32_t> SVnew;
            std::vector<uint32_t> SUnew;
            std::vector<uint32_t> SVPrime;
            std::unordered_set<uint32_t> SUSet;
            std::unordered_set<uint32_t> SVSet;
            for (uint32_t v : SV) {
                SVSet.insert(v);
                if (g->connectUV(pivot, v)) {
                    SVnew.push_back(v);

                } else {
                    SVPrime.push_back(v);
                }
            }

            for (uint32_t u : SU) {
                SUSet.insert(u);

                if (u != pivot) {
                    SUnew.push_back(u);
                }
            }
            SUnew.shrink_to_fit();
            SVnew.shrink_to_fit();
            NodeV* newNode = new NodeV(SUnew, SVnew, 1, 1, pivot);
            newNode->up += node->up + 1;
            newNode->vh = node->vh;
            newNode->uh = node->uh;
            newNode->vp = node->vp;
            // node->children.push_back(newNode);
            processNode(newNode);

            std::unordered_set<uint32_t> visited;

            for (uint32_t v_i : SVPrime) {
                visited.insert(v_i);
                std::vector<uint32_t> SVchild;
                std::vector<uint32_t> SUchild;
                std::vector<uint32_t> twoHopsCount(g->n2 + 1);

                for (uint32_t j = g->pV[v_i]; j < g->pV[v_i + 1]; j++) {
                    uint32_t u = g->e2[j];

                    if (SUSet.count(u) > 0) {
                        SUchild.push_back(u);
                    }

                    for (uint32_t k = g->pU[u]; k < g->pU[u + 1]; k++) {
                        uint32_t v = g->e1[k];

                        if (SVSet.count(v) > 0 && visited.count(v) == 0 && twoHopsCount[v] == 0) {
                            SVchild.push_back(v);
                        }
                        twoHopsCount[v]++;
                    }
                }

                NodeV* childNode = new NodeV(SUchild, SVchild, 2, 2, v_i);
                childNode->vh += node->vh + 1;
                childNode->uh = node->uh;
                childNode->vp = node->vp;
                childNode->up = node->up;
                // node->children.push_back(childNode);
                processNode(childNode);
            }
            // for (uint32_t v : SV) {
            //     SVBitMask[v] = 0;
            // }
            // for (uint32_t u : SU) {
            //     SUBitMask[u] = 0;
            // }

        } else {
            std::vector<uint32_t> SVnew;
            std::vector<uint32_t> SUnew;
            std::vector<uint32_t> SUPrime;
            std::unordered_set<uint32_t> SUSet;
            std::unordered_set<uint32_t> SVSet;
            for (uint32_t u : SU) {
                SUSet.insert(u);
                if (g->connectUV(u, pivot)) {
                    SUnew.push_back(u);
                } else {
                    SUPrime.push_back(u);
                }
            }

            for (uint32_t v : SV) {
                SVSet.insert(v);
                if (v != pivot) {
                    SVnew.push_back(v);
                }
            }

            NodeV* newNode = new NodeV(SUnew, SVnew, 1, 2, pivot);

            newNode->vp += node->vp + 1;
            newNode->vh = node->vh;
            newNode->uh = node->uh;
            newNode->up = node->up;
            processNode(newNode);

            std::unordered_set<uint32_t> visited;

            for (uint32_t u_i : SUPrime) {
                visited.insert(u_i);
                std::vector<uint32_t> SUchild;
                std::vector<uint32_t> SVchild;
                std::vector<uint32_t> twoHopsCount(g->n1 + 1);

                for (uint32_t j = g->pU[u_i]; j < g->pU[u_i + 1]; j++) {
                    uint32_t v = g->e1[j];
                    if (SVSet.count(v) > 0) {
                        SVchild.push_back(v);
                    }
                    for (uint32_t k = g->pV[v]; k < g->pV[v + 1]; k++) {
                        uint32_t u = g->e2[k];
                        if (SUSet.count(u) > 0 && visited.count(u) == 0 && twoHopsCount[u] == 0) {
                            SUchild.push_back(u);
                        }
                        twoHopsCount[u]++;
                    }
                }

                NodeV* childNode = new NodeV(SUchild, SVchild, 2, 1, u_i);
                // node->children.push_back(childNode);
                childNode->uh += node->uh + 1;
                childNode->vp = node->vp;
                childNode->vh = node->vh;
                childNode->up = node->up;
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

    if (SU.empty() || SV.empty()) {
        std::vector<uint32_t> SVnew;
        std::vector<uint32_t> SUnew;
        int side = 1;
        int vertex = 0;
        int pCount = 0;
        if (!SU.empty()) {
            pCount = SU.size();

        } else {
            pCount = SV.size();
            side = 2;
        }
        NodeV* newNode = new NodeV(SUnew, SVnew, 1, side, vertex, pCount);
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
        // node->children.push_back(newNode);
        processNode(newNode);
    }
}

std::pair<uint32_t, bool> BCTV2::selectPivot(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV) {
    auto Start = std::chrono::high_resolution_clock::now();
    uint32_t bestVertex = 0;
    bool isFromSU = true;
    uint32_t maxSum = 0;

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
    auto End = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(End - Start);

    return {bestVertex, isFromSU};
}