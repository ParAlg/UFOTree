#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "types.h"
#include "util.h"


template<typename aug_t>
struct UFOCluster {
    // UFO cluster data
    std::unordered_map<UFOCluster<aug_t>*,aug_t> neighbors;
    aug_t value;            // Stores subtree values or cluster path values
    UFOCluster<aug_t>* parent;
    // Constructor
    UFOCluster(aug_t value) : neighbors(), parent(), value(value) {};
    // Helper functions
    int get_degree();
    bool high_degree();
    bool parent_high_fanout();
    bool contracts();
    bool contains_neighbor(UFOCluster<aug_t>* c);
    void insert_neighbor(UFOCluster<aug_t>* c, aug_t value);
    void remove_neighbor(UFOCluster<aug_t>* c);
    UFOCluster<aug_t>* get_root();
};

template<typename aug_t>
class UFOTree {
public:
    // UFO tree interface
    UFOTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t dval);
    void link(vertex_t u, vertex_t v, aug_t value = 1);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
private:
    // Class data and parameters
    parlay::sequence<UFOCluster<aug_t>> leaves;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    parlay::sequence<std::unordered_set<UFOCluster<aug_t>*>> root_clusters;
    std::vector<std::pair<std::pair<UFOCluster<aug_t>*,UFOCluster<aug_t>*>,UFOCluster<aug_t>*>> contractions;
    // Helper functions
    void remove_ancestors(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v);
    void recluster_tree();
};

template<typename aug_t>
UFOTree<aug_t>::UFOTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) :
query_type(q), f(f), identity(id), default_value(d) {
    leaves.resize(n, d);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(6);
}

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename aug_t>
void UFOTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    remove_ancestors(&leaves[u], &leaves[v]);
    leaves[u].insert_neighbor(&leaves[v], value);
    leaves[v].insert_neighbor(&leaves[u], value);
    recluster_tree();
}

/* Cut vertex u and vertex v in the tree. */
template<typename aug_t>
void UFOTree<aug_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    remove_ancestors(&leaves[u], &leaves[v]);
    leaves[u].remove_neighbor(&leaves[v]);
    leaves[v].remove_neighbor(&leaves[u]);
    recluster_tree();
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename aug_t>
bool UFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

template<typename aug_t>
void UFOTree<aug_t>::remove_ancestors(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v) {

    int level = 0;
    auto prev_u = u;
    auto prev_v = v;
    auto curr_u = u->parent;
    auto curr_v = v->parent;
    bool del_u = false;
    bool del_v = false;

    while (curr_u && curr_v) {
        if (curr_u == curr_v) break;
        // U
        if (!curr_u->high_degree() && !prev_u->parent_high_fanout()) {  // We will delete curr_u next round
            for (auto entry : prev_u->neighbors) {
                auto neighbor = entry.first;
                if (neighbor->parent == prev_u->parent) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[0].insert(neighbor); // Keep track of parentless cluster
                }
            }
            if (del_u) { // Possibly delete prev_u
                for (auto entry : prev_u->neighbors)
                    entry.first->remove_neighbor(prev_u); // Remove prev from adjacency
                delete prev_u;
                root_clusters[level].erase(prev_u);
            }
            else prev_u->parent = nullptr;
            del_u = true;
        } else {                                                        // We will not delete curr_u next round
            if (del_u) { // Possibly delete prev_u
                for (auto entry : prev_u->neighbors)
                    entry.first->remove_neighbor(prev_u); // Remove prev from adjacency
                delete prev_u;
                root_clusters[level].erase(prev_u);
            }
            del_u = false;
        }
        // V
        if (!curr_v->high_degree() && !prev_v->parent_high_fanout()) {  // We will delete curr_v next round
            for (auto entry : prev_v->neighbors) {
                auto neighbor = entry.first;
                if (neighbor->parent == prev_v->parent) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[0].insert(neighbor); // Keep track of parentless cluster
                }
            }
            if (del_v) { // Possibly delete prev_v
                for (auto entry : prev_v->neighbors)
                    entry.first->remove_neighbor(prev_v); // Remove prev from adjacency
                delete prev_v;
                root_clusters[level].erase(prev_v);
            }
            else prev_v->parent = nullptr;
            del_v = true;
        } else {                                                        // We will not delete curr_v next round
            if (del_v) { // Possibly delete prev_v
                for (auto entry : prev_v->neighbors)
                    entry.first->remove_neighbor(prev_v); // Remove prev from adjacency
                delete prev_v;
                root_clusters[level].erase(prev_v);
            }
            del_v = false;
        }
        // Update pointers
        prev_u = curr_u;
        curr_u = prev_u->parent;
        prev_v = curr_v;
        curr_v = prev_v->parent;
    }
    // DO LAST DELETIONS
    if (curr_u == curr_v) curr_v = nullptr;
    if (!curr_u) {
        if (del_u) { // Possibly delete prev_u
            for (auto entry : prev_u->neighbors)
                entry.first->remove_neighbor(prev_u); // Remove prev from adjacency
            delete prev_u;
            root_clusters[level].erase(prev_u);
        }
    }
    if (!curr_v) {
        if (del_v) { // Possibly delete prev_v
            for (auto entry : prev_v->neighbors)
                entry.first->remove_neighbor(prev_v); // Remove prev from adjacency
            delete prev_v;
            root_clusters[level].erase(prev_v);
        }
    }
    // LOOP FOR REMAINING ONE PATH
    auto curr = curr_u;
    auto prev = prev_u;
    bool del = del_u;
    while (curr) {
        if (!curr->high_degree() && !prev->parent_high_fanout()) {  // We will delete curr next round
            for (auto entry : prev->neighbors) {
                auto neighbor = entry.first;
                if (neighbor->parent == prev->parent) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[0].insert(neighbor); // Keep track of parentless cluster
                }
            }
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry.first->remove_neighbor(prev); // Remove prev from adjacency
                delete prev;
                root_clusters[level].erase(prev);
            }
            else prev->parent = nullptr;
            del = true;
        } else {                                                        // We will not delete curr next round
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry.first->remove_neighbor(prev); // Remove prev from adjacency
                delete prev;
                root_clusters[level].erase(prev);
            }
            del = false;
        }
        prev = curr;
        curr = prev->parent;
    }
    // DO LAST DELETION
    if (!curr) {
        if (del) { // Possibly delete prev
            for (auto entry : prev->neighbors)
                entry.first->remove_neighbor(prev); // Remove prev from adjacency
            delete prev;
            root_clusters[level].erase(prev);
        }
    }
}

template<typename aug_t>
void UFOTree<aug_t>::recluster_tree() {
    return;
}

template<typename aug_t>
int UFOCluster<aug_t>::get_degree() {
    return neighbors.size();
}

template<typename aug_t>
bool UFOCluster<aug_t>::high_degree() {
    return (neighbors.size() > 2);
}

template<typename aug_t>
bool UFOCluster<aug_t>::parent_high_fanout() {
    assert(parent);
    if (neighbors.size() == 1)
        if (neighbors.begin()->first->get_degree() - parent->get_degree() > 1) return true;
    else if (this->get_degree() - parent->get_degree() > 1) return true;
    return false;
}

template<typename aug_t>
bool UFOCluster<aug_t>::contracts() {
    bool contracts = false;
    for (auto neighbor : this->neighbors)
        if (neighbor.first->parent == this->parent)
            contracts = true;
    return contracts;
}

template<typename aug_t>
bool UFOCluster<aug_t>::contains_neighbor(UFOCluster<aug_t>* c) {
    return (neighbors.find(c) != neighbors.end());
}

template<typename aug_t>
void UFOCluster<aug_t>::insert_neighbor(UFOCluster<aug_t>* c, aug_t value) {
    neighbors.insert({c, value});
}

template<typename aug_t>
void UFOCluster<aug_t>::remove_neighbor(UFOCluster<aug_t>* c) {
    neighbors.erase(c);
}

template<typename aug_t>
UFOCluster<aug_t>* UFOCluster<aug_t>::get_root() {
    UFOCluster<aug_t>* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}
