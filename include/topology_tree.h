#include <parlay/sequence.h>
#include "types.h"


template<typename aug_t>
struct TopologyCluster {
    // Topology cluster data
    TopologyCluster* neighbors[3];
    TopologyCluster* parent;
    aug_t value;
    // Helper functions
    void insert_neighbor(TopologyCluster<aug_t>* c);
    void remove_neighbor(TopologyCluster<aug_t>* c);
};

template<typename aug_t>
class TopologyTree {
public:
    // Topology Tree interface
    TopologyTree(vertex_t n, QueryType query_type, std::function<aug_t(aug_t, aug_t)> f);
    void link(vertex_t u, vertex_t v, aug_t value = 0);
    void cut(vertex_t u, vertex_t v);
private:
    // Class data and parameters
    parlay::sequence<TopologyCluster<aug_t>> leaves;
    std::function<aug_t(aug_t, aug_t)> f;
    QueryType query_type;
    parlay::sequence<parlay::sequence<TopologyCluster<aug_t>*>> root_clusters;
    // Helper functions
    void remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f) :
query_type(query_type), f(f) { leaves.reserve(n); }

template<typename aug_t>
void TopologyTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    remove_ancestors(u, v);
    leaves[u].insert_neighbor(&leaves[v]);
    leaves[v].insert_neighbor(&leaves[u]);
    // recluster(u, v);
}

template<typename aug_t>
void TopologyTree<aug_t>::cut(vertex_t u, vertex_t v) {
    remove_ancestors(u, v);
    leaves[u].remove_neighbor(&leaves[v]);
    leaves[v].remove_neighbor(&leaves[u]);
    // recluster(u, v);
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    // Collect all the root clusters while removing ancestors
    root_clusters.clear();
    root_clusters.reserve(HEIGHT(leaves.size()));
    for (auto neighbor : u->neighbors)
        if (neighbor->parent == u->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].push_back(neighbor); // Keep track of parentless cluster
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
            if (neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of parentless cluster
            }
            neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
    }
    for (auto neighbor : v->neighbors)
        if (neighbor->parent == v->parent) {
            neighbor->parent = nullptr; // Set sibling parent pointer to null
            root_clusters[0].push_back(neighbor); // Keep track of parentless cluster
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
            if (neighbor->parent == prev->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of parentless cluster
            }
            neighbor->remove_neighbor(prev); // Remove prev from adjacency
        }
        free(prev); // Remove cluster prev
    }
}

template<typename aug_t>
void TopologyCluster<aug_t>::insert_neighbor(TopologyCluster<aug_t>* c) {
    for (int i = 0; i < 3; i++)
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            return;
        }
    std::cerr << "No space to insert neighbor." << std::endl;
    std::abort();
}

template<typename aug_t>
void TopologyCluster<aug_t>::remove_neighbor(TopologyCluster<aug_t>* c) {
    for (int i = 0; i < 3; i++)
        if (this->neighbors[i] == c) {
            this->neighbors[i] = nullptr;
            return;
        }
    std::cerr << "Neighbor to delete not found." << std::endl;
    std::abort();
}
