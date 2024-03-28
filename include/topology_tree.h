#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "types.h"
#include "util.h"


template<typename aug_t>
struct TopologyCluster {
    // Topology cluster data
    TopologyCluster<aug_t>* neighbors[3];
    aug_t edge_values[3];   // Only for path queries
    aug_t value;            // Stores subtree values or cluster path values
    TopologyCluster<aug_t>* parent;
    // Constructor
    TopologyCluster(aug_t value) : neighbors(), edge_values(), parent(), value(value) {};
    // Helper functions
    int get_degree();
    bool contracts();
    bool contains_neighbor(TopologyCluster<aug_t>* c);
    void insert_neighbor(TopologyCluster<aug_t>* c, aug_t value);
    void remove_neighbor(TopologyCluster<aug_t>* c);
    TopologyCluster<aug_t>* get_root();
};

template<typename aug_t>
class TopologyTree {
public:
    // Topology tree interface
    TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t dval);
    void link(vertex_t u, vertex_t v, aug_t value = 1);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    aug_t subtree_query(vertex_t v, vertex_t p = MAX_VERTEX_T);
    aug_t path_query(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
private:
    // Class data and parameters
    parlay::sequence<TopologyCluster<aug_t>> leaves;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    parlay::sequence<std::unordered_set<TopologyCluster<aug_t>*>> root_clusters;
    std::vector<std::pair<std::pair<TopologyCluster<aug_t>*,TopologyCluster<aug_t>*>,TopologyCluster<aug_t>*>> contractions;
    // Helper functions
    void remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v);
    void recluster_tree();
};

template<typename aug_t>
TopologyTree<aug_t>::TopologyTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) :
query_type(q), f(f), identity(id), default_value(d) {
    leaves.resize(n, d);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(4);
}

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename aug_t>
void TopologyTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    remove_ancestors(&leaves[u], &leaves[v]);
    leaves[u].insert_neighbor(&leaves[v], value);
    leaves[v].insert_neighbor(&leaves[u], value);
    recluster_tree();
}

/* Cut vertex u and vertex v in the tree. */
template<typename aug_t>
void TopologyTree<aug_t>::cut(vertex_t u, vertex_t v) {
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
bool TopologyTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

/* Returns the value of the associative function f applied over
the augmented values for all the vertices in the subtree rooted
at v with respect to its parent p. If p = -1 (MAX_VERTEX_T) then
return the sum over the entire tree containing v. */
template<typename aug_t>
aug_t TopologyTree<aug_t>::subtree_query(vertex_t v, vertex_t p) {
    assert(v >= 0 && v < leaves.size() && p >= 0 && (p < leaves.size() || p == MAX_VERTEX_T));
    if (p == MAX_VERTEX_T) return leaves[v].get_root()->value;
    assert(leaves[v].contains_neighbor(&leaves[p]));
    // Get the total up until the LCA of v and p
    aug_t total = f(identity, leaves[v].value);
    auto curr_v = &leaves[v];
    auto curr_p = &leaves[p];
    while (curr_v->parent != curr_p->parent) {
        for (auto neighbor : curr_v->neighbors) if (neighbor && neighbor->parent == curr_v->parent)
            total = f(total, neighbor->value); // Only count vertices on the side of v
        curr_v = curr_v->parent;
        curr_p = curr_p->parent;
    }
    // Add the total after the LCA of v and p
    if (curr_v->get_degree() == 1) return total;
    if (curr_v->get_degree() == 3) {
        curr_v = curr_v->parent;
        while (curr_v->parent) {
            for (auto neighbor : curr_v->neighbors) if (neighbor && neighbor->parent == curr_v->parent)
                total = f(total, neighbor->value); // Count all remaining root clusters
            curr_v = curr_v->parent;
        }
        return total;
    }
    TopologyCluster<aug_t>* curr_u;
    for (auto neighbor : curr_v->neighbors) if (neighbor && neighbor != curr_p)
        curr_u = neighbor; // Find the neighbor of curr_v that is not curr_p
    while (curr_v->parent) {
        if (curr_u->parent == curr_v->parent) {
            total = f(total, curr_u->value); // Count the remaining root clusters on the side of v
            if (curr_u->get_degree() == 1) return total;
            if (curr_u->get_degree() == 3) {
                curr_v = curr_v->parent;
                while (curr_v->parent) {
                    for (auto neighbor : curr_v->neighbors) if (neighbor && neighbor->parent == curr_v->parent)
                        total = f(total, neighbor->value); // Count all remaining root clusters
                    curr_v = curr_v->parent;
                }
                return total;
            }
            for (auto neighbor : curr_u->neighbors) if (neighbor && neighbor != curr_v)
                curr_u = neighbor->parent; // Find the neighbor of curr_u that is not curr_v
            curr_v = curr_v->parent;
        } else {
            curr_u = curr_u->parent;
            curr_v = curr_v->parent;
        }
    }
    return total;
}

/* Returns the value of the associative function f applied over
the augmented values for all the edges on the unique path from
vertex u to vertex v. */
template<typename aug_t>
aug_t TopologyTree<aug_t>::path_query(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && connected(u,v));
    // Compute the path on both sides for both vertices until they combine
    aug_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = path_v1 = path_v2 = identity;
    TopologyCluster<aug_t> *bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
    bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = nullptr;
    if (leaves[u].get_degree() == 2) {
        bdry_u1 = leaves[u].neighbors[0] ? leaves[u].neighbors[0] : leaves[u].neighbors[1];
        bdry_u2 = leaves[u].neighbors[2] ? leaves[u].neighbors[2] : leaves[u].neighbors[1];
    }
    if (leaves[v].get_degree() == 2) {
        bdry_v1 = leaves[v].neighbors[0] ? leaves[v].neighbors[0] : leaves[v].neighbors[1];
        bdry_v2 = leaves[v].neighbors[2] ? leaves[v].neighbors[2] : leaves[v].neighbors[1];
    }
    auto curr_u = &leaves[u];
    auto curr_v = &leaves[v];
    while (curr_u->parent != curr_v->parent) {
        for (int i = 0; i < 3; i++) {
            auto neighbor = curr_u->neighbors[i];
            if (neighbor && neighbor->parent == curr_u->parent) {
                if (curr_u->get_degree() == 2) {
                    if (curr_u->parent->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_u1) {
                            path_u1 = f(path_u1, f(curr_u->edge_values[i], neighbor->value));
                            bdry_u2 = bdry_u2->parent;
                            for (int i = 0; i < 3; i++)
                                if(curr_u->parent->neighbors[i] && curr_u->parent->neighbors[i] != bdry_u2)
                                    bdry_u1 = curr_u->parent->neighbors[i];
                        } else {
                            path_u2 = f(path_u2, f(curr_u->edge_values[i], neighbor->value));
                            bdry_u1 = bdry_u1->parent;
                            for (int i = 0; i < 3; i++)
                                if(curr_u->parent->neighbors[i] && curr_u->parent->neighbors[i] != bdry_u1)
                                    bdry_u2 = curr_u->parent->neighbors[i];
                        }
                    } else {
                        // Binary to Unary
                        path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
                    }
                } else {
                    if (curr_u->parent->get_degree() == 2) {
                        // Unary to Binary
                        path_u1 = path_u2 = f(path_u1, curr_u->edge_values[i]);
                        bdry_u1 = curr_u->parent->neighbors[0] ? curr_u->parent->neighbors[0] : curr_u->parent->neighbors[1];
                        bdry_u2 = curr_u->parent->neighbors[2] ? curr_u->parent->neighbors[2] : curr_u->parent->neighbors[1];
                    } else {
                        // Unary to Unary
                        path_u1 = f(path_u1, f(curr_u->edge_values[i], neighbor->value));
                    }
                }
                break;
            }
        }
        if (!curr_u->contracts()) {
            if (bdry_u1) bdry_u1 = bdry_u1->parent;
            if (bdry_u2) bdry_u2 = bdry_u2->parent;
        }
        curr_u = curr_u->parent;
        for (int i = 0; i < 3; i++) {
            auto neighbor = curr_v->neighbors[i];
            if (neighbor && neighbor->parent == curr_v->parent) {
                if (curr_v->get_degree() == 2) {
                    if (curr_v->parent->get_degree() == 2) {
                        // Binary to Binary
                        if (neighbor == bdry_v1) {
                            path_v1 = f(path_v1, f(curr_v->edge_values[i], neighbor->value));
                            bdry_v2 = bdry_v2->parent;
                            for (int i = 0; i < 3; i++)
                                if(curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v2)
                                    bdry_v1 = curr_v->parent->neighbors[i];
                        } else {
                            path_v2 = f(path_v2, f(curr_v->edge_values[i], neighbor->value));
                            bdry_v1 = bdry_v1->parent;
                            for (int i = 0; i < 3; i++)
                                if(curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v1)
                                    bdry_v2 = curr_v->parent->neighbors[i];
                        }
                    } else {
                        // Binary to Unary
                        path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                    }
                } else {
                    if (curr_v->parent->get_degree() == 2) {
                        // Unary to Binary
                        path_v1 = path_v2 = f(path_v1, curr_v->edge_values[i]);
                        bdry_v1 = curr_v->parent->neighbors[0] ? curr_v->parent->neighbors[0] : curr_v->parent->neighbors[1];
                        bdry_v2 = curr_v->parent->neighbors[2] ? curr_v->parent->neighbors[2] : curr_v->parent->neighbors[1];
                    } else {
                        // Unary to Unary
                        path_v1 = f(path_v1, f(curr_v->edge_values[i], neighbor->value));
                    }
                }
                break;
            }
        }
        if (!curr_v->contracts()) {
            if (bdry_v1) bdry_v1 = bdry_v1->parent;
            if (bdry_v2) bdry_v2 = bdry_v2->parent;
        }
        curr_v = curr_v->parent;
    }
    // Get the correct path sides when the two vertices meet at the LCA
    aug_t total = identity;
    if (curr_u->get_degree() == 2)
        total = f(total, (curr_v == bdry_u1) ? path_u1 : path_u2);
    else
        total = f(total, path_u1);
    if (curr_v->get_degree() == 2)
        total = f(total, (curr_u == bdry_v1) ? path_v1 : path_v2);
    else
        total = f(total, path_v1);
    // Add the value of the last edge
    for (int i = 0; i < 3; i++) if (curr_u->neighbors[i] == curr_v)
        total = f(total, curr_u->edge_values[i]);
    return total;
}

template<typename aug_t>
void TopologyTree<aug_t>::remove_ancestors(TopologyCluster<aug_t>* u, TopologyCluster<aug_t>* v) {
    for (auto leaf: {u,v}) {
        for (auto neighbor : leaf->neighbors) {
            if (neighbor && neighbor->parent == leaf->parent) {
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
                            parent = new TopologyCluster<aug_t>(default_value);
                            root_clusters[level+1].insert(parent);
                        }
                        if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                for (int i = 0; i < 3; i++) {
                    auto neighbor = cluster->neighbors[i];
                    if (neighbor && !neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        auto parent = new TopologyCluster<aug_t>(default_value);
                        if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
                        if (query_type == PATH && neighbor->get_degree() == 2)
                            parent->value = f(f(cluster->value, neighbor->value), cluster->edge_values[i]);
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
                for (int i = 0; i < 3; i++) {
                    auto neighbor = cluster->neighbors[i];
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        if (neighbor->contracts()) continue;
                        auto parent = neighbor->parent;
                        if (!parent) parent = new TopologyCluster<aug_t>(default_value);
                        if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
                        if (query_type == PATH && neighbor->get_degree() == 2)
                            parent->value = f(f(cluster->value, neighbor->value), cluster->edge_values[i]);
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
                        auto parent = new TopologyCluster<aug_t>(default_value);
                        if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                        if (!parent) parent = new TopologyCluster<aug_t>(default_value);
                        if (query_type == SUBTREE) parent->value = f(cluster->value, neighbor->value);
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
                auto parent = new TopologyCluster<aug_t>(default_value);
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
                    parent->insert_neighbor(c1->neighbors[i]->parent, c1->edge_values[i]);
                    c1->neighbors[i]->parent->insert_neighbor(parent, c1->edge_values[i]);
                }
            }
            for (int i = 0; i < 3; ++i) {
                if (c2->neighbors[i] && c2->neighbors[i] != c1) { // Don't add c1's parent (itself)
                    parent->insert_neighbor(c2->neighbors[i]->parent, c2->edge_values[i]);
                    c2->neighbors[i]->parent->insert_neighbor(parent, c2->edge_values[i]);
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
void TopologyCluster<aug_t>::insert_neighbor(TopologyCluster<aug_t>* c, aug_t value) {
    if (this->contains_neighbor(c)) return;
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            this->edge_values[i] = value;
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
