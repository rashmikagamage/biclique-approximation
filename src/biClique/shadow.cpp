#include <vector>
#include <set>
#include <iostream>

struct ShadowPair {
    std::set<uint32_t> SU;
    std::set<uint32_t> SV;
};

class shadow {
public:
    void ShadowConstruction() {

        std::vector<ShadowPair> shadow;

        // Iterate over all vertices in U
        for (uint32_t u = 0; u < g->n1; u++) {
            // Iterate over all edges from vertex u
            for (uint32_t i = g->pU[u]; i < g->pU[u + 1]; i++) {
                uint32_t v = g->e1[i]; // v is a neighbor of u

                // Get neighbors of u excluding v
                std::set<uint32_t> neighborsU;
                for (uint32_t j = g->pU[u]; j < g->pU[u + 1]; j++) {
                    if (g->e1[j] != v) {
                        neighborsU.insert(g->e1[j]);
                    }
                }

                // Get neighbors of v excluding u
                std::set<uint32_t> neighborsV;
                for (uint32_t j = g->pV[v]; j < g->pV[v + 1]; j++) {
                    if (g->e2[j] != u) {
                        neighborsV.insert(g->e2[j]);
                    }
                }

                // Find common neighbors of neighborsV
                std::set<uint32_t> commonNeighbors;
                for (uint32_t w : neighborsV) {
                    if (neighborsU.find(w) != neighborsU.end()) {
                        commonNeighbors.insert(w);
                    }
                }

                // Create a shadow pair and add to the shadow list
                ShadowPair sp;
                sp.SU = neighborsU;
                sp.SV = commonNeighbors;
                shadow.push_back(sp);
            }
        }

        // Output or process the shadow pairs
        for (const auto& sp : shadow) {
            std::cout << "Shadow SU: ";
            for (uint32_t u : sp.SU) {
                std::cout << u << " ";
            }
            std::cout << "\nShadow SV: ";
            for (uint32_t v : sp.SV) {
                std::cout << v << " ";
            }
            std::cout << std::endl;
        }
    }
};