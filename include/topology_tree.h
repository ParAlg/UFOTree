#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "types.h"
#include "util.h"


template<typename aug_t>
struct TopologyCluster {
    // Topology cluster data
    TopologyCluster<aug_t>* neighbors[3];
    TopologyCluster<aug_t>* parent;
    aug_t value;
    // Helper functions
    int get_degree();
    bool contains_neighbor(TopologyCluster<aug_t>* c);
    void insert_neighbor(TopologyCluster<aug_t>* c);
    void remove_neighbor(TopologyCluster<aug_t>* c);
    TopologyCluster<aug_t>* get_root();
};

template<typename aug_t>
class TopologyTree {
public:
    // Topology Tree interface
    TopologyTree(vertex_t n, QueryType query_type, std::function<aug_t(aug_t, aug_t)> f);
    void link(vertex_t u, vertex_t v, aug_t value = 0);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
private:
    // Class data and parameters
    parlay::sequence<TopologyCluster<aug_t>> leaves;
    std::function<aug_t(aug_t, aug_t)> f;
    QueryType query_type;
    parlay::sequence<parlay::sequence<TopologyCluster<aug_t>*>> root_clusters;
    // Helper functions
    void remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
    void recluster_tree();
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f) :
query_type(query_type), f(f) { leaves.resize(n); }

template<typename aug_t>
void TopologyTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    remove_ancestors(&leaves[u], &leaves[v]);
    leaves[u].insert_neighbor(&leaves[v]);
    leaves[v].insert_neighbor(&leaves[u]);
    recluster_tree();
}

template<typename aug_t>
void TopologyTree<aug_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    remove_ancestors(&leaves[u], &leaves[v]);
    leaves[u].remove_neighbor(&leaves[v]);
    leaves[v].remove_neighbor(&leaves[u]);
    recluster_tree();
}

template<typename aug_t>
bool TopologyTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    // Collect all the root clusters while removing ancestors
    root_clusters.clear();
    root_clusters.reserve(max_tree_height(leaves.size()));
    for (auto neighbor : u->neighbors) {
        if (neighbor && neighbor->parent == u->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].push_back(neighbor); // Keep track of parentless cluster
        }
    }
    auto curr = u->parent;
    int level = 0;
    u->parent = nullptr;
    root_clusters[0].push_back(u);
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        level++;
        for (auto neighbor : prev->neighbors) {
            if (neighbor && neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of parentless cluster
            }
            if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
    }
    for (auto neighbor : v->neighbors) {
        if (neighbor && neighbor->parent == v->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].push_back(neighbor); // Keep track of parentless cluster
        }
    }
    curr = v->parent;
    level = 0;
    v->parent = nullptr;
    root_clusters[0].push_back(v);
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        level++;
        for (auto neighbor : prev->neighbors) {
            if (neighbor && neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of parentless cluster
            }
            if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
    }
    level = -1;
    while (root_clusters[++level].size() > 0) root_clusters[level] = parlay::remove_duplicates(root_clusters[level]);
}

template<typename aug_t>
void TopologyTree<aug_t>::recluster_tree() {
    parlay::sequence<std::pair<std::pair<TopologyCluster<aug_t>*,TopologyCluster<aug_t>*>,TopologyCluster<aug_t>*>> contractions;
    int level = -1;
    while (root_clusters[++level].size() > 0) {
        contractions.clear();
        // Combine deg 3 root clusters with deg 1 root clusters
        for (auto cluster : root_clusters[level]) {
            if (cluster->get_degree() == 3) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 1 or 2 root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->parent == nullptr && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 2 non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 2 && neighbor->parent != nullptr) {
                        bool neighbor_contracts = false; // Check if neighbor already contracts
                        for (auto grandneighbor : neighbor->neighbors)
                            if (grandneighbor && grandneighbor->parent == neighbor->parent)
                                neighbor_contracts = true;
                        if (neighbor_contracts) continue;
                        auto parent = neighbor->parent;
                        if (!parent) parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Combine deg 1 root clusters with deg 1 root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Combine deg 1 root clusters with deg 2 or 3 non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->parent != nullptr && (neighbor->get_degree() == 2 || neighbor->get_degree() == 3)) {
                        bool neighbor_contracts = false; // Check if neighbor already contracts
                        for (auto grandneighbor : neighbor->neighbors)
                            if (grandneighbor && grandneighbor->parent == neighbor->parent)
                                neighbor_contracts = true;
                        if (neighbor_contracts) continue;
                        auto parent = neighbor->parent;
                        if (!parent) parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Add remaining uncombined clusters to the next level
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() > 0) {
                auto parent = new TopologyCluster<aug_t>();
                parent->value = cluster->value;
                cluster->parent = parent;
                root_clusters[level+1].push_back(parent);
                contractions.push_back({{cluster,cluster},parent});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1 = contraction.first.first;
            auto c2 = contraction.first.second;
            auto parent = contraction.second;
            for (int i = 0; i < 3; i++) parent->neighbors[i] = nullptr;
            for (int i = 0; i < 3; i++) {
                if (c1->neighbors[i] && c1->neighbors[i] != c2) { // Don't add c2's parent (itself)
                    parent->insert_neighbor(c1->neighbors[i]->parent);
                    c1->neighbors[i]->parent->insert_neighbor(parent);
                }
            }
            for (int i = 0; i < 3; i++) if (c2->neighbors[i]) {
                if (c2->neighbors[i] && c2->neighbors[i] != c1) { // Don't add c1's parent (itself)
                    parent->insert_neighbor(c2->neighbors[i]->parent);
                    c2->neighbors[i]->parent->insert_neighbor(parent);
                }
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
    }
}

template<typename aug_t>
int TopologyCluster<aug_t>::get_degree() {
    int deg = 0;
    for (auto neighbor : this->neighbors) if (neighbor) deg++;
    return deg;
}


template<typename aug_t>
bool TopologyCluster<aug_t>::contains_neighbor(TopologyCluster<aug_t>* c) {
    for (int i = 0; i < 3; i++) if (this->neighbors[i] == c) return true;
    return false;
}

template<typename aug_t>
void TopologyCluster<aug_t>::insert_neighbor(TopologyCluster<aug_t>* c) {
    if (this->contains_neighbor(c)) return;
    for (int i = 0; i < 3; i++) {
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            return;
        }
    }
    std::cerr << "No space to insert neighbor." << std::endl;
    std::abort();
}

template<typename aug_t>
void TopologyCluster<aug_t>::remove_neighbor(TopologyCluster<aug_t>* c) {
    for (int i = 0; i < 3; i++) {
        if (this->neighbors[i] == c) {
            this->neighbors[i] = nullptr;
            return;
        }
    }
    std::cerr << "Neighbor to delete not found." << std::endl;
    std::abort();
}

template<typename aug_t>
TopologyCluster<aug_t>* TopologyCluster<aug_t>::get_root() {
    TopologyCluster<aug_t>* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}
