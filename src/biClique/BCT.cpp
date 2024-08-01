#include "BCT.h"

#include <signal.h>

#include <algorithm>
#include <iostream>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>
int current_i = 0;

void segfault_handler(int signal) {
    printf("Caught segfault at %d\n", current_i);
    exit(1);
}

void BCT::buildTree() {
    signal(SIGSEGV, segfault_handler);
    Node root;

    std::vector<uint32_t> SU(g->n1 + 1);
    std::vector<uint32_t> SV(g->n2 + 1);
    std::stack<Node*> Q;
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
        if (SU.size() == 0 && SV.size() == 0) {
            continue;
        }

        Node* n = new Node(SU, SV, 2, 1, u);
        root.children.push_back(n);
        Q.push(n);
    }

    while (!Q.empty()) {
        Node* node = Q.top();

        Q.pop();
        SU.clear();
        SV.clear();
        SU = node->SU;
        SV = node->SV;
        // printf("-------------------------------------------------\n");
        // printf("Q size: %d\n", Q.size());
        // printf("node->SU size: %d\n", node->SU.size());
        // printf("node->SV size: %d\n", node->SV.size());
        // printf("SU: ");

        // for (int i = 0; i < SU.size(); i++) {
        //     printf("%d ", SU[i]);
        // }

        // printf("\n");
        // printf("SV: ");
        // for (int i = 0; i < SV.size(); i++) {
        //     printf("%d ", SV[i]);
        // }
        // printf("\n");
        if (SU.size() == 0 && SV.size() == 0) {
            continue;
        }

        if (SV.size() != 0 && SU.size() != 0) {
            std::pair<uint32_t, bool> pivotPair = selectPivot(SU, SV);
            // std::pair<uint32_t, bool> pivotPair = selectPivoteWithSide(SU, SV, node->vertexSide);
            ui pivot = pivotPair.first;

            // pivot from U
            if (pivotPair.second) {
                std::vector<uint32_t> SVnew;
                std::vector<uint32_t> SUnew;

                std::vector<uint32_t> SVPrime;

                for (uint32_t v : SV) {
                    if (g->connectUV(pivot, v)) {
                        SVnew.push_back(v);
                    } else {
                        SVPrime.push_back(v);
                    }
                }

                // Calculate S_u_new as SU \ pivot
                for (uint32_t u : SU) {
                    if (u != pivot) {
                        SUnew.push_back(u);
                    }
                }

                // printf("pivot: %d\n", pivot);
                // printf("SUnew: ");
                // for (int i = 0; i < SUnew.size(); i++) {
                //     printf("%d ", SUnew[i]);
                // }
                // printf("\n");
                // printf("SVnew: ");
                // for (int i = 0; i < SVnew.size(); i++) {
                //     printf("%d ", SVnew[i]);
                // }
                // printf("\n");
                // printf("SVPrime: ");
                // for (int i = 0; i < SVPrime.size(); i++) {
                //     printf("%d ", SVPrime[i]);
                // }
                // printf("\n");

                // printf("node.vertex: %d\n", node->vertex);
                // printf("node.vertexSide: %d\n", node->vertexSide);
                // printf("node.label: %d\n", node->label);

                Node* newNode = new Node(SUnew, SVnew, 1, 1, pivot);
                node->children.push_back(newNode);
                Q.push(newNode);

                // biGraph sg = createSubgraph(SUnew, SV);
                std::unordered_set<uint32_t> visited;

                for (uint32_t i = 0; i < SVPrime.size(); i++) {
                    uint32_t v_i = SVPrime[i];
                    visited.insert(v_i);
                    std::vector<uint32_t> SVchild;
                    std::vector<uint32_t> SUchild;
                    std::vector<uint32_t> twoHopsCount(g->n2 + 1);

                    for (uint32_t j = g->pV[v_i]; j < g->pV[v_i + 1]; j++) {
                        uint32_t u = g->e2[j];
                        // auto it = std::find(SU.begin(), SU.end(), u);
                        // if (it != SU.end()) {
                        //     SUchild.push_back(u);
                        // }
                        current_i = u;
                        for (uint32_t k = g->pU[u]; k < g->pU[u + 1]; k++) {
                            uint32_t v = g->e1[k];
                            // auto it = std::find(SV.begin(), SV.end(), v);
                            // if (it != SV.end()) {
                            //     if (visited.count(v) == 0 && twoHopsCount[v] == 0) {
                            //         SVchild.push_back(v);
                            //     }
                            // }
                            twoHopsCount[v]++;
                        }
                    }

                    Node* childNode = new Node(SUchild, SVchild, 2, 2, v_i);
                    node->children.push_back(childNode);
                    Q.push(childNode);
                    // printf("SUchild: ");
                    // for (int i = 0; i < SUchild.size(); i++) {
                    //     printf("%d ", SUchild[i]);
                    // }
                    // printf("\n");
                    // printf("SVchild: ");
                    // for (int i = 0; i < SVchild.size(); i++) {
                    //     printf("%d ", SVchild[i]);
                    // }
                    // printf("\n");
                }
            } else {
                // If pivot is from SV
                printf("Pivot is from SV\n");
                fflush(stdout);
                std::vector<uint32_t> SVnew;
                std::vector<uint32_t> SUnew;

                // Calculate S_u_new as S_u intersect N(rho)
                std::vector<uint32_t> SUPrime;

                for (uint32_t u : SU) {
                    if (g->connectUV(u, pivot)) {
                        SUnew.push_back(u);
                    } else {
                        SUPrime.push_back(u);
                    }
                }

                // Calculate S_v_new as SV \ rho
                for (uint32_t v : SV) {
                    if (v != pivot) {
                        SVnew.push_back(v);
                    }
                }

                Node* newNode = new Node(SU, SVnew, 1, 2);
                node->children.push_back(newNode);
                Q.push(newNode);

                biGraph sg = createSubgraph(SUPrime, SVnew);

                std::unordered_set<uint32_t> visited;
                for (ui i = 0; i < SUPrime.size(); i++) {
                    uint32_t u_i = SUPrime[i];
                    visited.insert(u_i);
                    std::vector<uint32_t> SUchild;
                    std::vector<uint32_t> SVchild;
                    std::vector<uint32_t> twoHopCount(g->n1 + 1);
                    ui u_iMap = sg.SUMap[u_i];
                    for (ui j = sg.pU[u_iMap]; j < sg.pU[u_iMap + 1]; j++) {
                        uint32_t v = sg.e1[j];
                        SVchild.push_back(v);
                        ui vMap = sg.SVMap[v];
                        for (ui k = sg.pV[vMap]; k < sg.pV[vMap + 1]; k--) {
                            uint32_t u = sg.e2[k];
                            if (twoHopCount[u] == 0 && visited.find(u) == visited.end()) {
                                SUchild.push_back(u);
                            }
                            twoHopCount[u]++;
                        }
                    }

                    Node* childNode = new Node(SUchild, SVchild, 2, 2);
                    node->children.push_back(childNode);
                    Q.push(childNode);
                }
            }
        }
        if (SU.size() == 0 || SV.size() == 0) {
            std::vector<uint32_t> SVnew;
            std::vector<uint32_t> SUnew;
            int side = 0;
            int vertex = 0;
            if (SU.size() == 0) {
                side = 2;
                vertex = SV[0];
            } else {
                side = 1;
                vertex = SU[0];
            }

            for (int i = 1; i < SU.size(); i++) {
                SUnew.push_back(SU[i]);
            }
            for (int i = 1; i < SV.size(); i++) {
                SVnew.push_back(SV[i]);
            }
            // printf("one one one\n");
            // printf("SUnew: ");
            // for (int i = 0; i < SUnew.size(); i++) {
            //     printf("%d ", SUnew[i]);
            // }
            // printf("\n");
            // printf("SVnew: ");
            // for (int i = 0; i < SVnew.size(); i++) {
            //     printf("%d ", SVnew[i]);
            // }
            // printf("\n");
            Node* newNode = new Node(SUnew, SVnew, 1, side, vertex);
            node->children.push_back(newNode);
            Q.push(newNode);
        }
    }
    printf("Tree building completed\n");

    std::vector<std::vector<double>> result = countBicliques(&root);
    for (int i = 2; i < 6; i++) {
        for (int j = 2; j < 6; j++) {
            printf("(%d - %d) = %f\n", i, j, result[i][j]);
        }
    }
}

std::pair<uint32_t, bool> BCT::selectPivot(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV) {
    uint32_t bestVertex = 0;
    bool isFromSU = true;
    uint32_t maxSum = 0;

    for (uint32_t u : SU) {
        std::vector<uint32_t> PU, PV;
        for (uint32_t v : SV) {
            if (g->connectUV(u, v)) {
                PV.push_back(v);
            }
        }

        uint32_t sum = SU.size() - 1 + PV.size();
        if (sum >= maxSum) {
            maxSum = sum;
            bestVertex = u;
            isFromSU = true;
        }
    }

    // for (uint32_t v : SV) {
    //     std::vector<uint32_t> PU, PV;
    //     for (uint32_t u : SU) {
    //         if (g->connectUV(u, v)) {
    //             PU.push_back(u);
    //         }
    //     }
    //     uint32_t sum = PU.size() + SV.size() - 1;
    //     if (sum >= maxSum) {
    //         maxSum = sum;
    //         bestVertex = v;
    //         isFromSU = false;
    //     }
    // }

    return {SU[SU.size() - 1], true};
}

std::pair<uint32_t, bool> BCT::selectPivoteWithSide(std::vector<uint32_t>& SU, std::vector<uint32_t>& SV, int side) {
    uint32_t bestVertex = 0;
    bool isFromSU = true;
    uint32_t maxSum = 0;

    if (side == 2) {
        for (uint32_t u : SU) {
            std::vector<uint32_t> PU, PV;
            for (uint32_t v : SV) {
                if (g->connectUV(u, v)) {
                    PV.push_back(v);
                }
            }

            uint32_t sum = SU.size() - 1 + PV.size();
            if (sum >= maxSum) {
                maxSum = sum;
                bestVertex = u;
                isFromSU = true;
            }
        }
    } else {
        for (uint32_t v : SV) {
            std::vector<uint32_t> PU, PV;
            for (uint32_t u : SU) {
                if (g->connectUV(u, v)) {
                    PU.push_back(u);
                }
            }
            uint32_t sum = PU.size() + SV.size() - 1;
            if (sum >= maxSum) {
                maxSum = sum;
                bestVertex = v;
                isFromSU = false;
            }
        }
    }

    return {bestVertex, isFromSU};
}
biGraph BCT::createSubgraph(const std::vector<uint32_t>& minSU, const std::vector<uint32_t>& minSV) {
    biGraph sg;

    uint32_t llSize = minSU.size();
    uint32_t rrSize = minSV.size();

    sg.n1 = llSize;
    sg.n2 = rrSize;
    sg.pU.resize(llSize + 2, 0);
    sg.pV.resize(rrSize + 2, 0);
    sg.e1.resize(llSize * rrSize + 1);
    sg.e2.resize(llSize * rrSize + 1);
    for (uint32_t i = 0; i < llSize; i++) {
        sg.SUMap[minSU[i]] = i;
    }
    for (uint32_t i = 0; i < rrSize; i++) {
        sg.SVMap[minSV[i]] = i;
    }

    int edgeIndex = 0;
    for (uint32_t j = 0; j < llSize; j++) {
        uint32_t x = minSU[j];
        sg.pU[j + 1] = sg.pU[j];
        for (uint32_t k = 0; k < rrSize; k++) {
            uint32_t y = minSV[k];
            if (g->connectUV(x, y)) {
                sg.e1[edgeIndex] = y;
                edgeIndex++;
                sg.pU[j + 1]++;
                sg.pV[k + 1]++;
            }
        }
    }

    sg.e1.resize(edgeIndex);

    for (uint32_t v = 0; v < rrSize; v++) {
        sg.pV[v + 1] += sg.pV[v];
    }

    edgeIndex = 0;

    for (uint32_t u = 0; u < llSize; u++) {
        for (uint32_t i = sg.pU[u]; i < sg.pU[u + 1]; i++) {
            uint32_t v = sg.e1[i];

            int index = sg.SVMap[v];
            sg.e2[sg.pV[index]] = minSU[u];

            sg.pV[index]++;

            edgeIndex++;
        }
    }
    for (int v = rrSize; v >= 1; v--) {
        sg.pV[v] = sg.pV[v - 1];
    }
    assert(sg.PV[minSV.size()] == sg.pU[minSU.size()]);
    sg.pV[0] = 0;

    return sg;
}

void BCT::traversePaths(Node* node, std::vector<Node*>& path, std::vector<std::vector<Node*>>& allPaths) {
    path.push_back(node);
    if (node->children.empty()) {
        allPaths.push_back(path);
    } else {
        for (Node* child : node->children) {
            traversePaths(child, path, allPaths);
        }
    }

    path.pop_back();
}

std::vector<std::vector<double>> BCT::countBicliques(Node* root) {
    int buf = g->maxDu * g->maxDv + 1;
    std::vector<std::vector<double>> count(buf, std::vector<double>(buf, 0));
    std::vector<std::vector<Node*>> allPaths;
    std::vector<Node*> path;

    // Traverse all root to leaf paths
    traversePaths(root, path, allPaths);

    for (const std::vector<Node*>& path : allPaths) {
        int up = 0, uh = 0, vp = 0, vh = 0;

        for (Node* node : path) {
            if (node->vertexSide == 1) {
                if (node->label == 1) {
                    up++;
                } else if (node->label == 2) {
                    uh++;
                }
            } else if (node->vertexSide == 2) {
                if (node->label == 1) {
                    vp++;
                } else if (node->label == 2) {
                    vh++;
                }
            }
        }
        for (int i = 0; i <= up; ++i) {
            for (int j = 0; j <= vp; ++j) {
                current_i = j + vh;
                count[i + uh][j + vh] += C[up][i] * C[vp][j];
                // count[i][j] += C[up + uh][i] * C[vp + vh][j];
            }
        }
    }

    return count;
}