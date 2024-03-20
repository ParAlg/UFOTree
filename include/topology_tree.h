#include <unordered_set>
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
    bool contracts();
    bool contains_neighbor(TopologyCluster<aug_t>* c);
    void insert_neighbor(TopologyCluster<aug_t>* c);
    void remove_neighbor(TopologyCluster<aug_t>* c);
    TopologyCluster<aug_t>* get_root();
};

template<typename aug_t>
class TopologyTree {
public:
    // Topology tree interface
    TopologyTree(vertex_t n, QueryType query_type, std::function<aug_t(aug_t, aug_t)> f);
    void link(vertex_t u, vertex_t v, aug_t value = 0);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
private:
    // Class data and parameters
    parlay::sequence<TopologyCluster<aug_t>> leaves;
    std::function<aug_t(aug_t, aug_t)> f;
    QueryType query_type;
    parlay::sequence<std::unordered_set<TopologyCluster<aug_t>*>> root_clusters;
    std::vector<std::pair<std::pair<TopologyCluster<aug_t>*,TopologyCluster<aug_t>*>,TopologyCluster<aug_t>*>> contractions;
    // Helper functions
    void remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
    void recluster_tree();
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f) :
query_type(query_type), f(f) { leaves.resize(n); root_clusters.resize(max_tree_height(n)); contractions.reserve(4); }

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
    for (auto neighbor : u->neighbors) {
        if (neighbor && neighbor->parent == u->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].insert(neighbor); // Keep track of parentless cluster
        }
    }
    auto curr = u->parent;
    int level = 0;
    u->parent = nullptr;
    root_clusters[0].insert(u);
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        level++;
        for (auto neighbor : prev->neighbors) {
            if (neighbor && neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].insert(neighbor); // Keep track of parentless cluster
            }
            if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
        root_clusters[level].erase(prev);
    }
    for (auto neighbor : v->neighbors) {
        if (neighbor && neighbor->parent == v->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].insert(neighbor); // Keep track of parentless cluster
        }
    }
    curr = v->parent;
    level = 0;
    v->parent = nullptr;
    root_clusters[0].insert(v);
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        level++;
        for (auto neighbor : prev->neighbors) {
            if (neighbor && neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].insert(neighbor); // Keep track of parentless cluster
            }
            if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
        root_clusters[level].erase(prev);
    }
}

template<typename aug_t>
void TopologyTree<aug_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
        if (root_clusters[level].empty()) continue;
        // Combine deg 3 root clusters with deg 1 root  or non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (cluster->get_degree() == 3) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto parent = neighbor->parent;
                        if (!parent) { // If neighbor is a root cluster
                            parent = new TopologyCluster<aug_t>();
                            root_clusters[level+1].insert(parent);
                        }
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
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
                    if (neighbor && !neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent = new TopologyCluster<aug_t>();
                        parent->value = f(cluster->value, neighbor->value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].insert(parent);
                        contractions.push_back({{cluster,neighbor},parent});
                        break;
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        if (neighbor->contracts()) continue;
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
                        root_clusters[level+1].insert(parent);
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
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 2 || neighbor->get_degree() == 3)) {
                        if (neighbor->contracts()) continue;
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
                root_clusters[level+1].insert(parent);
                contractions.push_back({{cluster,cluster},parent});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1 = contraction.first.first;
            auto c2 = contraction.first.second;
            auto parent = contraction.second;
            for (int i = 0; i < 3; ++i) parent->neighbors[i] = nullptr;
            for (int i = 0; i < 3; ++i) {
                if (c1->neighbors[i] && c1->neighbors[i] != c2) { // Don't add c2's parent (itself)
                    parent->insert_neighbor(c1->neighbors[i]->parent);
                    c1->neighbors[i]->parent->insert_neighbor(parent);
                }
            }
            for (int i = 0; i < 3; ++i) {
                if (c2->neighbors[i] && c2->neighbors[i] != c1) { // Don't add c1's parent (itself)
                    parent->insert_neighbor(c2->neighbors[i]->parent);
                    c2->neighbors[i]->parent->insert_neighbor(parent);
                }
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename aug_t>
int TopologyCluster<aug_t>::get_degree() {
    int deg = 0;
    for (auto neighbor : this->neighbors) if (neighbor) deg++;
    return deg;
}

template<typename aug_t>
bool TopologyCluster<aug_t>::contracts() {
    bool contracts = false;
    for (auto neighbor : this->neighbors)
        if (neighbor && neighbor->parent == this->parent)
            contracts = true;
    return contracts;
}

template<typename aug_t>
bool TopologyCluster<aug_t>::contains_neighbor(TopologyCluster<aug_t>* c) {
    for (int i = 0; i < 3; ++i) if (this->neighbors[i] == c) return true;
    return false;
}

template<typename aug_t>
void TopologyCluster<aug_t>::insert_neighbor(TopologyCluster<aug_t>* c) {
    if (this->contains_neighbor(c)) return;
    for (int i = 0; i < 3; ++i) {
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
    for (int i = 0; i < 3; ++i) {
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
