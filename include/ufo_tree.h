#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "types.h"
#include "util.h"


template<typename aug_t>
struct UFOCluster {
    // UFO cluster data
    std::unordered_map<UFOCluster<aug_t>*,aug_t> neighbors;
    aug_t value; // Stores subtree values or cluster path values
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
    void print_tree();
private:
    // Class data and parameters
    parlay::sequence<UFOCluster<aug_t>> leaves;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    parlay::sequence<std::unordered_set<UFOCluster<aug_t>*>> root_clusters;
    std::vector<std::pair<std::pair<UFOCluster<aug_t>*,UFOCluster<aug_t>*>,bool>> contractions;
    // Helper functions
    void remove_ancestors(UFOCluster<aug_t>* c);
    void recluster_tree();
    void disconnect_siblings(UFOCluster<aug_t>* c, int level);
    void insert_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v, aug_t value);
    void remove_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v);
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
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    insert_adjacency(&leaves[u], &leaves[v], value);
    recluster_tree();
}

/* Cut vertex u and vertex v in the tree. */
template<typename aug_t>
void UFOTree<aug_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    remove_adjacency(&leaves[u], &leaves[v]);
    recluster_tree();
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename aug_t>
bool UFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

/* Removes the ancestors of cluster c that are not high degree nor
high fan-out and add them to root_clusters. */
template<typename aug_t>
void UFOTree<aug_t>::remove_ancestors(UFOCluster<aug_t>* c) {
    int level = 0;
    auto prev = c;
    auto curr = c->parent;
    bool del = false;
    while (curr) {
        if (!curr->high_degree() && !prev->parent_high_fanout()) { // We will delete curr next round
            disconnect_siblings(prev, level);
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry.first->remove_neighbor(prev); // Remove prev from adjacency
                delete prev;
                root_clusters[level].erase(prev);
            } else {
                prev->parent = nullptr;
                root_clusters[level].insert(prev);
            }
            del = true;
        } else { // We will not delete curr next round
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry.first->remove_neighbor(prev); // Remove prev from adjacency
                delete prev;
                root_clusters[level].erase(prev);
            } else if (prev->get_degree() <= 1) {
                prev->parent = nullptr;
                root_clusters[level].insert(prev);
            }
            del = false;
        }
        // Update pointers
        prev = curr;
        curr = prev->parent;
        level++;
    }
    // DO LAST DELETIONS
    if (!curr) {
        if (del) { // Possibly delete prev
            for (auto entry : prev->neighbors)
                entry.first->remove_neighbor(prev); // Remove prev from adjacency
            delete prev;
            root_clusters[level].erase(prev);
        } else root_clusters[level].insert(prev);
    }
}

template<typename aug_t>
void UFOTree<aug_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
        if (root_clusters[level].empty()) continue;
        // Merge deg exactly 3 root clusters with all of its deg 1 neighbors
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 3) {
                auto parent = new UFOCluster<aug_t>(default_value);
                cluster->parent = parent;
                root_clusters[level+1].insert(parent);
                bool first_contraction = true;
                for (auto entry : cluster->neighbors)
                    if (entry.first->get_degree() == 1) {
                        auto old_parent = entry.first->parent;
                        entry.first->parent = cluster->parent;
                        while (old_parent) {
                            auto temp = old_parent;
                            old_parent = old_parent->parent;
                            delete temp;
                        }
                        contractions.push_back({{entry.first, cluster},first_contraction});
                        first_contraction = false;
                    }
                if (first_contraction) contractions.push_back({{cluster, cluster}, true});
            }
        }
        // Always combine deg 1 root clusters with its neighboring cluster
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) {
                auto neighbor = cluster->neighbors.begin()->first;  // Deg 1 cluster only one neighbor
                if (neighbor->parent) {
                    if (neighbor->get_degree() == 2 && neighbor->contracts()) continue;
                    cluster->parent = neighbor->parent;
                    contractions.push_back({{cluster,neighbor},false});
                } else {
                    auto parent = new UFOCluster<aug_t>(default_value);
                    cluster->parent = parent;
                    neighbor->parent = parent;
                    root_clusters[level+1].insert(parent);
                    contractions.push_back({{cluster,neighbor},true});
                }
                if (neighbor->get_degree() == 3) { // For deg exactly 3 make all deg 1 neighbors combine with it
                    for (auto entry : neighbor->neighbors) {
                        if (entry.first->get_degree() == 1 && entry.first->parent != neighbor->parent) {
                            auto old_parent = entry.first->parent;
                            entry.first->parent = neighbor->parent;
                            contractions.push_back({{entry.first, neighbor},false});
                            neighbor->parent->remove_neighbor(old_parent);
                            if (old_parent) delete old_parent;
                        }
                    }
                }
                if (cluster->parent->get_degree() == 1 || cluster->parent->get_degree() == 2) {
                    auto prev = cluster->parent;
                    auto curr = cluster->parent->parent;
                    bool del = (prev->parent && !prev->contracts());
                    while (del) {
                        del = (curr->parent && !curr->contracts());
                        for (auto entry : curr->neighbors) entry.first->remove_neighbor(curr);
                        delete curr;
                        prev = curr;
                        curr = prev->parent;
                    }
                    cluster->parent->parent = nullptr;
                    root_clusters[level+1].insert(cluster->parent);
                }
            }
        }
        // Combine deg 2 root clusters with deg 2 root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto entry : cluster->neighbors) {
                    auto neighbor = entry.first;
                    if (!neighbor->parent && (neighbor->get_degree() == 2)) {
                        auto parent = new UFOCluster<aug_t>(default_value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].insert(parent);
                        contractions.push_back({{cluster,neighbor},true});
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
                        cluster->parent = neighbor->parent;
                        contractions.push_back({{cluster,neighbor},false}); // The order here is important
                        auto prev = cluster->parent;
                        auto curr = cluster->parent->parent;
                        bool del = (prev->parent && !prev->contracts());
                        if (!del) break;
                        while (del) {
                            del = (curr->parent && !curr->contracts());
                            for (auto entry : curr->neighbors) entry.first->remove_neighbor(curr);
                            delete curr;
                            prev = curr;
                            curr = prev->parent;
                        }
                        cluster->parent->parent = nullptr;
                        root_clusters[level+1].insert(cluster->parent);
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
                contractions.push_back({{cluster,cluster},true});
            }
        }
        // Fill in the neighbor lists of the new clusters
        for (auto contraction : contractions) {
            auto c1 = contraction.first.first;
            auto c2 = contraction.first.second;
            auto parent = c1->parent;
            bool new_parent = contraction.second;
            if (new_parent) {
                for (auto entry : c1->neighbors) {
                    if (entry.first != c2) {
                        parent->insert_neighbor(entry.first->parent, entry.second);
                        entry.first->parent->insert_neighbor(parent, entry.second);
                    }
                }
                if (c1 != c2)
                for (auto entry : c2->neighbors) {
                    if (entry.first != c1) {
                        parent->insert_neighbor(entry.first->parent, entry.second);
                        entry.first->parent->insert_neighbor(parent, entry.second);
                    }
                }
            } else {
                if (c1->get_degree() == 1) parent->remove_neighbor(c1->parent);
                if (c2->get_degree() == 1) parent->remove_neighbor(c2->parent);
                if (c1->get_degree() == 2) // We ordered contractions so c2 is the one that had a parent already
                    for (auto entry : c1->neighbors) if (entry.first != c2)
                        insert_adjacency(parent, entry.first->parent, entry.second);
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

template<typename aug_t>
void UFOTree<aug_t>::disconnect_siblings(UFOCluster<aug_t>* c, int level) {
    if (c->get_degree() == 1) {
        auto center = c->neighbors.begin()->first;
        if (center->parent != c->parent) return;
        for (auto entry : center->neighbors) {
            auto neighbor = entry.first;
            if (neighbor->parent == c->parent && neighbor != c) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].insert(neighbor); // Keep track of root clusters
            }
        }
        center->parent = nullptr;
        root_clusters[level].insert(center);
    } else {
        for (auto entry : c->neighbors) {
            auto neighbor = entry.first;
            if (neighbor->parent == c->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].insert(neighbor); // Keep track of root clusters
            }
        }
    }
}

template<typename aug_t>
void UFOTree<aug_t>::insert_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v, aug_t value) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->insert_neighbor(curr_v, value);
        curr_v->insert_neighbor(curr_u, value);
        curr_u = curr_u->parent;
        curr_v = curr_v->parent;
    }
}

template<typename aug_t>
void UFOTree<aug_t>::remove_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->remove_neighbor(curr_v);
        curr_v->remove_neighbor(curr_u);
        curr_u = curr_u->parent;
        curr_v = curr_v->parent;
    }
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
    // if (neighbors.size() == 1)
    //     if (neighbors.begin()->first->get_degree() - parent->get_degree() > 1) return true;
    // else if (this->get_degree() - parent->get_degree() > 1) return true;
    // return false;

    if (get_degree() == 1) {
        if (neighbors.begin()->first->get_degree() - parent->get_degree() > 2) return true;
    } else {
        if (get_degree() - parent->get_degree() > 2) return true;
    }
    return false;
}

template<typename aug_t>
bool UFOCluster<aug_t>::contracts() {
    for (auto neighbor : this->neighbors)
        if (neighbor.first->parent == this->parent)
            return true;
    return false;
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
