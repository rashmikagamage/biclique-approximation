#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <queue>
#include <unordered_set>
#include <vector>

struct BiGraph {
    uint32_t n1, n2;
    std::vector<uint32_t> deg1, deg2, pU, pV, e1, e2;

    bool connectUV(uint32_t u, uint32_t v) const {
        return std::find(e1.begin() + pU[u], e1.begin() + pU[u + 1], v) != e1.begin() + pU[u + 1];
    }
};

struct Node {
    std::vector<uint32_t> Sv, Su;
    std::string edgeLabel;
    std::vector<Node*> children;

    Node(const std::vector<uint32_t>& Sv, const std::vector<uint32_t>& Su, const std::string& edgeLabel)
        : Sv(Sv), Su(Su), edgeLabel(edgeLabel) {}
};

void printTree(const Node* node, const std::string& prefix = "", bool isLast = true) {
    std::cout << prefix << (isLast ? "└─" : "├─") << " Node: Su={";
    for (const auto& u : node->Su) {
        std::cout << u << " ";
    }
    std::cout << "} Sv={";
    for (const auto& v : node->Sv) {
        std::cout << v << " ";
    }
    std::cout << "} EdgeLabel=" << node->edgeLabel << "\n";

    for (size_t i = 0; i < node->children.size(); ++i) {
        printTree(node->children[i], prefix + (isLast ? "   " : "│  "), i == node->children.size() - 1);
    }
}

Node* buildBCTree(BiGraph* g) {
    Node* root = new Node(std::vector<uint32_t>(g->n2), std::vector<uint32_t>(g->n1), "root");
    std::iota(root->Su.begin(), root->Su.end(), 0);
    std::iota(root->Sv.begin(), root->Sv.end(), 0);

    std::queue<Node*> Q;
    Q.push(root);

    while (!Q.empty()) {
        Node* current = Q.front();
        Q.pop();

        if (current->Sv.empty() && current->Su.empty()) {
            continue;
        }

        if (!current->Sv.empty() && !current->Su.empty()) {
            // Select pivot
            uint32_t rho = current->Su[0];
            for (uint32_t u : current->Su) {
                if (g->deg1[u] > g->deg1[rho]) {
                    rho = u;
                }
            }

            std::vector<uint32_t> S_v_prime, S_u_prime;
            for (uint32_t j = g->pU[rho]; j < g->pU[rho + 1]; ++j) {
                uint32_t v = g->e1[j];
                if (std::find(current->Sv.begin(), current->Sv.end(), v) != current->Sv.end()) {
                    S_v_prime.push_back(v);
                }
            }

            for (uint32_t u : current->Su) {
                if (u != rho) {
                    S_u_prime.push_back(u);
                }
            }

            Node* child1 = new Node(S_v_prime, S_u_prime, "(" + std::to_string(rho) + ", p)");
            current->children.push_back(child1);
            Q.push(child1);

            std::vector<uint32_t> remainingSv;
            for (uint32_t v : current->Sv) {
                if (v != rho) {
                    remainingSv.push_back(v);
                }
            }

            for (uint32_t v : remainingSv) {
                std::vector<uint32_t> newSv, newSu;

                for (uint32_t j = g->pV[v]; j < g->pV[v + 1]; ++j) {
                    uint32_t u = g->e2[j];
                    newSu.push_back(u);
                }

                for (uint32_t u : newSu) {
                    for (uint32_t k = g->pU[u]; k < g->pU[u + 1]; ++k) {
                        uint32_t w = g->e1[k];
                        if (std::find(newSv.begin(), newSv.end(), w) == newSv.end()) {
                            newSv.push_back(w);
                        }
                    }
                }

                Node* child2 = new Node(newSv, newSu, "(" + std::to_string(v) + ", h)");
                current->children.push_back(child2);
                Q.push(child2);
            }
        } else {
            if (!current->Sv.empty()) {
                for (uint32_t v : current->Sv) {
                    auto S_v_prime = current->Sv;
                    S_v_prime.erase(std::remove(S_v_prime.begin(), S_v_prime.end(), v), S_v_prime.end());
                    Node* child = new Node(S_v_prime, std::vector<uint32_t>(), "(" + std::to_string(v) + ", p)");
                    current->children.push_back(child);
                    Q.push(child);
                }
            }

            if (!current->Su.empty()) {
                for (uint32_t u : current->Su) {
                    auto S_u_prime = current->Su;
                    S_u_prime.erase(std::remove(S_u_prime.begin(), S_u_prime.end(), u), S_u_prime.end());
                    Node* child = new Node(std::vector<uint32_t>(), S_u_prime, "(" + std::to_string(u) + ", p)");
                    current->children.push_back(child);
                    Q.push(child);
                }
            }
        }
    }

    return root;
}

int main() {
    std::ifstream infile("../data/dense.txt");
    if (!infile) {
        std::cerr << "Unable to open file input.txt";
        return 1;
    }

    uint32_t n, m, e;
    infile >> n >> m >> e;

    BiGraph g;
    g.n1 = n;
    g.n2 = m;
    g.deg1.resize(n, 0);
    g.deg2.resize(m, 0);
    g.pU.resize(n + 1, 0);
    g.pV.resize(m + 1, 0);

    std::vector<std::pair<uint32_t, uint32_t>> edges;
    for (uint32_t i = 0; i < e; ++i) {
        uint32_t u, v;
        infile >> u >> v;
        edges.emplace_back(u, v);
        g.deg1[u]++;
        g.deg2[v]++;
    }

    g.pU[0] = 0;
    for (uint32_t i = 1; i <= n; ++i) {
        g.pU[i] = g.pU[i - 1] + g.deg1[i - 1];
    }

    g.pV[0] = 0;
    for (uint32_t i = 1; i <= m; ++i) {
        g.pV[i] = g.pV[i - 1] + g.deg2[i - 1];
    }

    g.e1.resize(e);
    g.e2.resize(e);
    std::vector<uint32_t> indexU(n, 0), indexV(m, 0);
    for (const auto& edge : edges) {
        uint32_t u = edge.first, v = edge.second;
        g.e1[g.pU[u] + indexU[u]] = v;
        g.e2[g.pV[v] + indexV[v]] = u;
        indexU[u]++;
        indexV[v]++;
    }

    Node* treeRoot = buildBCTree(&g);

    // Print the tree
    printTree(treeRoot);

    return 0;
}
