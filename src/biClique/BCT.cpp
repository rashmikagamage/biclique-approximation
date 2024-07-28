#include "BCT.h"

#include <algorithm>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

void addChild(Node& parent, Node& child) {
    parent.children.push_back(child);
}

void BCT::buildTree() {
    Node root;
    std::vector<uint32_t> SU(g->n1 + 1);
    std::vector<uint32_t> SV(g->n2 + 1);
    std::queue<Node> Q;
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
        Node n(SU, SV, 2);
        createSubgraph(SU, SV);

        addChild(root, n);
        Q.push(n);
    }

    while (!Q.empty()) {
        Node node = Q.front();
        Q.pop();
        SU.clear();
        std::vector<uint32_t> SU = node.SU;
        std::vector<uint32_t> SV = node.SV;

        if (SU.size() == 0 && SV.size() == 0) {
            continue;
        }

        if (SV.size() != 0 && SU.size() != 0) {
            std::pair<uint32_t, bool> pivotPair = selectPivot(SU, SV);
            ui pivot = pivotPair.first;

            if (pivotPair.second) {
                std::vector<uint32_t> SVnew;
                std::vector<uint32_t> SUnew;

                // Calculate S_v_new as S_v intersect N(rho)
                std::vector<uint32_t> SVPrime;

                for (uint32_t v : SV) {
                    if (g->connectUV(pivot, v)) {
                        SVnew.push_back(v);
                    } else {
                        SVPrime.push_back(v);
                    }
                }

                // Calculate S_u_new as SU \ rho
                for (uint32_t u : SU) {
                    if (u != pivot) {
                        SUnew.push_back(u);
                    }
                }

                Node newNode(SUnew, SVnew, 1);
                addChild(node, newNode);
                Q.push(newNode);

                biGraph sg = createSubgraph(SU, SV);

                // Create child nodes for each v_i in S_v_diff
                for (ui i = 0; i < SVPrime.size(); i++) {
                    uint32_t v_i = SVPrime[i];
                    std::vector<uint32_t> SVchild;
                    std::vector<uint32_t> SUchild;
                    std::vector<uint32_t> twoHopCount(sg.n2 + 1);
                    ui v_iMap = sg.SVMap[v_i];
                    for (ui j = sg.pV[v_iMap]; j < sg.pV[v_iMap + 1]; j++) {
                        uint32_t u = sg.e2[j];
                        SUchild.push_back(u);
                        ui uMap = sg.SUMap[u];
                        for (ui k = sg.pU[uMap + 1] - 1; k >= sg.pU[uMap]; k--) {
                            uint32_t v = sg.e1[k];
                            if (v == v_i) {
                                break;
                            }
                            if (twoHopCount[v] == 0) {
                                SVchild.push_back(v);
                            }
                            twoHopCount[v]++;
                        }
                    }

                    Node childNode(SUchild, SVchild, 2);
                    addChild(newNode, childNode);
                    Q.push(childNode);
                }

            } else {
                // If pivot is from SV
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

                Node newNode(SUnew, SVnew, 1);
                addChild(node, newNode);
                Q.push(newNode);

                biGraph sg = createSubgraph(SU, SV);

                for (ui i = 0; i < SUPrime.size(); i++) {
                    uint32_t u_i = SUPrime[i];
                    std::vector<uint32_t> SUchild;
                    std::vector<uint32_t> SVchild;
                    std::vector<uint32_t> twoHopCount(sg.n1 + 1);
                    ui u_iMap = sg.SUMap[u_i];
                    for (ui j = sg.pU[u_iMap]; j < sg.pU[u_iMap + 1]; j++) {
                        uint32_t v = sg.e1[j];
                        SVchild.push_back(v);
                        ui vMap = sg.SVMap[v];
                        for (ui k = sg.pV[vMap]; k < sg.pV[vMap + 1]; k--) {
                            uint32_t u = sg.e2[k];
                            if (u == u_i) {
                                break;
                            }
                            if (twoHopCount[u] == 0) {
                                SUchild.push_back(u);
                            }
                            twoHopCount[u]++;
                        }
                    }
                    Node childNode(SUchild, SVchild, 2);
                    addChild(newNode, childNode);
                    Q.push(childNode);
                }
            }
        } else {
            std::vector<uint32_t> SVnew;
            std::vector<uint32_t> SUnew;

            for (int i = 1; i < SU.size(); i++) {
                SUnew.push_back(SU[i]);
            }
            for (int i = 1; i < SV.size(); i++) {
                SVnew.push_back(SV[i]);
            }
            Node newNode(SUnew, SVnew, 1);
            addChild(node, newNode);
            Q.push(newNode);
        }
    }

    traversePaths(root);
    while (root.children.size() > 0) {
        root.children.pop_back();
    }
}

void BCT::printPath(const std::vector<int>& path) {
    for (int label : path) {
        std::cout << label << " ";
    }
    std::cout << std::endl;
}

void BCT::dfs(Node& node, std::vector<int>& currentPath) {
    currentPath.push_back(node.label);

    if (node.children.empty()) {
        printPath(currentPath);
    } else {
        for (Node& child : node.children) {
            dfs(child, currentPath);
        }
    }

    currentPath.pop_back();
}

void BCT::traversePaths(Node& root) {
    std::vector<int> currentPath;
    dfs(root, currentPath);
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
        if (sum > maxSum) {
            maxSum = sum;
            bestVertex = u;
            isFromSU = true;
        }
    }

    for (uint32_t v : SV) {
        std::vector<uint32_t> PU, PV;
        for (uint32_t u : SU) {
            if (g->connectUV(u, v)) {
                PU.push_back(u);
            }
        }
        uint32_t sum = PU.size() + SV.size() - 1;
        if (sum > maxSum) {
            maxSum = sum;
            bestVertex = v;
            isFromSU = false;
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

    sg.pV[0] = 0;

    return sg;
}