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
    for (auto leaf: {u,v}) {
        for (auto entry : leaf->neighbors) {
            auto neighbor = entry.first;
            if (neighbor->parent == leaf->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[0].insert(neighbor); // Keep track of parentless cluster
            }
        }
        auto curr = leaf->parent;
        int level = 0;
        leaf->parent = nullptr;
        root_clusters[0].insert(leaf);
        while (curr) {
            auto prev = curr;
            curr = prev->parent;
            level++;
            for (auto entry : prev->neighbors) {
                auto neighbor = entry.first;
                if (neighbor->parent == prev->parent) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].insert(neighbor); // Keep track of parentless cluster
                }
                neighbor->remove_neighbor(prev); // Remove prev from adjacency
            }
            free(prev); // Remove cluster prev
            root_clusters[level].erase(prev);
        }
        for (auto entry : v->neighbors) {
            auto neighbor = entry.first;
            if (neighbor->parent == v->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[0].insert(neighbor); // Keep track of parentless cluster
            }
        }
    }
}

template<typename aug_t>
void UFOTree<aug_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
        if (root_clusters[level].empty()) continue;
        // Combine deg 3 root clusters with deg 1 root  or non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (cluster->get_degree() == 3) {
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (neighbor->get_degree() == 1) {
                        auto parent = neighbor->parent;
                        if (!parent) { // If neighbor is a root cluster
                            parent = new UFOCluster<aug_t>(default_value);
                            root_clusters[level+1].insert(parent);
                        }
                        // if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (!neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent = new UFOCluster<aug_t>(default_value);
                        // if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
                        // if (query_type == PATH && neighbor->get_degree() == 2)
                        //     parent->value = f(f(cluster->value, neighbor->value), cluster->edge_values[i]);
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
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        if (neighbor->contracts()) continue;
                        auto parent = neighbor->parent;
                        if (!parent) parent = new UFOCluster<aug_t>(default_value);
                        // if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
                        // if (query_type == PATH && neighbor->get_degree() == 2)
                        //     parent->value = f(f(cluster->value, neighbor->value), cluster->edge_values[i]);
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
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (neighbor->get_degree() == 1) {
                        auto parent = new UFOCluster<aug_t>(default_value);
                        // if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (neighbor->parent && (neighbor->get_degree() == 2 || neighbor->get_degree() == 3)) {
                        if (neighbor->contracts()) continue;
                        auto parent = neighbor->parent;
                        if (!parent) parent = new UFOCluster<aug_t>(default_value);
                        // if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                auto parent = new UFOCluster<aug_t>(default_value);
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
            for (auto entry : c1->neighbors) {
                if (entry.first != c2) {
                    parent->insert_neighbor(entry.first->parent, entry.second);
                    entry.first->parent->insert_neighbor(parent, entry.second);
                }
            }
            for (auto entry : c2->neighbors) {
                if (entry.first != c1) {
                    parent->insert_neighbor(entry.first->parent, entry.second);
                    entry.first->parent->insert_neighbor(parent, entry.second);
                }
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename aug_t>
int UFOCluster<aug_t>::get_degree() {
    return neighbors.size();
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
