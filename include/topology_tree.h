#include <parlay/sequence.h>
#include "types.h"


template<typename aug_t>
struct TopologyCluster {
    TopologyCluster* neighbors[3];
    TopologyCluster* parent;
    aug_t value;
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
    // Helper functions
    void remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
    void insert_neighbor(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
    void remove_neighbor(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f) :
query_type(query_type), f(f) { leaves.reserve(n); }

template<typename aug_t>
void TopologyTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    remove_ancestors(u, v);
    insert_neighbor(&leaves[u], &leaves[v]);
    insert_neighbor(&leaves[v], &leaves[u]);
    // recluster(u, v);
}

template<typename aug_t>
void TopologyTree<aug_t>::cut(vertex_t u, vertex_t v) {
    remove_ancestors(u, v);
    remove_neighbor(&leaves[u], &leaves[v]);
    remove_neighbor(&leaves[v], &leaves[u]);
    // recluster(u, v);
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    parlay::sequence<std::array<TopologyCluster<aug_t>*,4>> root_clusters;
    auto curr = u->parent;
    for (auto neighbor : u->neighbors) // Set sibling parent pointer to null
        if (neighbor->parent == u->parent)
            neighbor->parent = nullptr;
    u->parent = nullptr;
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        for (auto neighbor : prev->neighbors) // Set sibling parent pointer to null
            if (neighbor->parent == prev->parent)
                neighbor->parent = nullptr;
        free(prev);
    }
    curr = v->parent;
    for (auto neighbor : v->neighbors) // Set sibling parent pointer to null
        if (neighbor->parent == v->parent)
            neighbor->parent = nullptr;
    while (curr) {
        auto prev = curr;
        curr = prev->parent;
        for (auto neighbor : prev->neighbors) // Set sibling parent pointer to null
            if (neighbor->parent == prev->parent)
                neighbor->parent = nullptr;
        free(prev);
    }
}

template<typename aug_t>
void TopologyTree<aug_t>::insert_neighbor(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    for (int i = 0; i < 3; i++)
        if (leaves[u]->neighbors[i] == nullptr) {
            leaves[u]->neighbors[i] = leaves[v];
            return;
        }
    std::cerr << "No space to insert neighbor." << std::endl;
    std::abort();
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_neighbor(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    for (int i = 0; i < 3; i++)
        if (leaves[u]->neighbors[i] == leaves[v]) {
            leaves[u]->neighbors[i] = nullptr;
            return;
        }
    std::cerr << "Neighbor to delete not found." << std::endl;
    std::abort();
}
