#pragma once
#include "types.h"
#include "util.h"
#include "ufo_cluster.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>


template<typename v_t, typename e_t>
class UFOTree {
using Cluster = UFOCluster<v_t, e_t>;
public:
    // UFO tree interface
    UFOTree(
        vertex_t n, QueryType q = CONNECTIVITY,
        std::function<v_t(v_t, v_t)> f_v = [](v_t x, v_t y) -> v_t {return x;},
        std::function<e_t(e_t, e_t)> f_e = [](e_t x, e_t y) -> e_t {return x;});
    UFOTree(
        vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e);
    UFOTree(int n, QueryType q, std::function<v_t(v_t, v_t)> f, v_t id, v_t d_val);
    ~UFOTree();
    void link(vertex_t u, vertex_t v);
    void link(vertex_t u, vertex_t v, e_t value);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    e_t path_query(vertex_t u, vertex_t v);
    // Testing helpers
    size_t space();
    size_t count_nodes();
    size_t get_height();
    bool is_valid();
    void print_tree();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<Cluster*> del_clusters;
    std::vector<Cluster*> root_clusters;
    std::vector<Cluster*> next_del_clusters;
    std::vector<Cluster*> next_root_clusters;
    UpdateType update_type;
    e_t link_value;
    Cluster *update_u, *update_v;
    // Query functions and parameters
    QueryType query_type;
    std::function<v_t(v_t, v_t)> f_v;
    v_t identity_v;
    v_t default_v;
    std::function<e_t(e_t, e_t)> f_e;
    e_t identity_e;
    e_t default_e;
    // We preallocate UFO clusters and store unused clusters in free_clusters
    std::vector<Cluster*> free_clusters;
    Cluster* allocate_cluster();
    void free_cluster(Cluster* c);
    // Helper functions
    void add_del(Cluster * c);
    void add_next_del(Cluster * c);
    bool should_delete(Cluster* cluster);
    void disconnect_children(Cluster* c);
    void update_tree();
    void recluster_roots();
};

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e)
    : query_type(q), f_v(f_v), f_e(f_e) {
    leaves.resize(n);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(vertex_t n, QueryType q,
        std::function<v_t(v_t, v_t)> f_v, std::function<e_t(e_t, e_t)> f_e,
        v_t id_v, e_t id_e, v_t dval_v, e_t dval_e)
    : query_type(q), f_v(f_v), f_e(f_e), identity_v(id_v), identity_e(id_e),
     default_v(dval_v), default_e(dval_e) {
    leaves.resize(n, default_v);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::UFOTree(int n, QueryType q,
        std::function<v_t(v_t, v_t)> f, v_t id, v_t d_val)
    : query_type(q), f_v(f), identity_v(id), default_v(d_val) {
    if constexpr (std::is_same<v_t,e_t>::value) {
        f_e = f;
        identity_e = id;
        default_e = d_val;
    }
    leaves.resize(n, default_v);
    for (int i = 0; i < n; ++i)
        free_clusters.push_back(new Cluster());
}

template<typename v_t, typename e_t>
UFOTree<v_t, e_t>::~UFOTree() {
    // Clear all memory
    std::unordered_set<Cluster*> clusters;
    for (auto leaf : leaves) {
        auto curr = leaf.parent;
        while (curr) {
            clusters.insert(curr);
            curr = curr->parent;
        }
    }
    for (auto cluster : clusters) delete cluster;
    for (auto cluster : free_clusters) delete cluster;
}

template<typename v_t, typename e_t>
UFOCluster<v_t, e_t>* UFOTree<v_t, e_t>::allocate_cluster() {
    if (!free_clusters.empty()) {
        auto c = free_clusters.back();
        free_clusters.pop_back();
        return c;
    }
    return new Cluster();
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::free_cluster(UFOCluster<v_t, e_t>* c) {
    c->parent = nullptr;
    if (c->has_neighbor_set()) [[unlikely]] delete c->get_neighbor_set();
    for (int i = 0; i < UFO_NEIGHBOR_MAX; ++i)
        c->neighbors[i] = nullptr;
    c->degree = 0;
    if (c->has_child_set()) [[unlikely]] delete c->get_child_set();
    for (int i = 0; i < UFO_CHILD_MAX; ++i)
        c->children[i] = nullptr;
    c->fanout = 0;
    c->mark = false;
    free_clusters.push_back(c);
}

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::link(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    update_type = INSERT;
    update_u = &leaves[u];
    update_v = &leaves[v];
    if (!leaves[u].parent) root_clusters.push_back(&leaves[u]);
    else add_del(leaves[u].parent);
    if (!leaves[v].parent) root_clusters.push_back(&leaves[v]);
    else add_del(leaves[v].parent);
    update_tree();
}
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::link(vertex_t u, vertex_t v, e_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    update_type = INSERT;
    link_value = value;
    update_u = &leaves[u];
    update_v = &leaves[v];
    if (!leaves[u].parent) root_clusters.push_back(&leaves[u]);
    else add_del(leaves[u].parent);
    if (!leaves[v].parent) root_clusters.push_back(&leaves[v]);
    else add_del(leaves[v].parent);
    update_tree();
}

/* Cut vertex u and vertex v in the tree. */
template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    update_type = DELETE;
    update_u = &leaves[u];
    update_v = &leaves[v];
    if (!leaves[u].parent) root_clusters.push_back(&leaves[u]);
    else add_del(leaves[u].parent);
    if (!leaves[v].parent) root_clusters.push_back(&leaves[v]);
    else add_del(leaves[v].parent);
    update_tree();
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::add_del(Cluster * c) {
    if (c && !c->mark) {
        c->mark = true;
        del_clusters.push_back(c);
    }
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::add_next_del(Cluster * c) {
    if (c && !c->mark) {
        c->mark = true;
        next_del_clusters.push_back(c);
    }
}

template<typename v_t, typename e_t>
inline bool UFOTree<v_t, e_t>::should_delete(Cluster* cluster) {
    if (cluster->degree > 2 || cluster->fanout > 2) return false;
    // if (cluster->fanout == 2) {
    //     auto child1 = cluster->children[0];
    //     auto child2 = cluster->children[1];
    //     if (child1->degree + child2->degree <= 4) return false;
    // }
    // if (cluster->fanout == 1) {
    //     auto child = cluster->children[0];
    //     bool can_contract = false;
    //     FOR_ALL_NEIGHBORS(child, [&](Cluster* neighbor, e_t _) {
    //         if (child->degree + neighbor->degree <= 4) can_contract = true;
    //     });
    //     if (can_contract) return true;
    // }
    return true;
}

template<typename v_t, typename e_t>
inline void UFOTree<v_t, e_t>::disconnect_children(Cluster* c) {
    FOR_ALL_CHILDREN(c, [&](Cluster* child) {
        child->parent = nullptr;
        root_clusters.push_back(child);
    });
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::update_tree() {
    while (!root_clusters.empty() || !del_clusters.empty()) {
        // For deletions delete the edge in level i
        if (update_type == DELETE && update_u && update_v) {
            update_u->remove_neighbor(update_v);
            update_v->remove_neighbor(update_u);
            if (update_u->parent == update_v->parent) {
                if (update_u->degree > update_v->degree)
                    std::swap(update_u, update_v);
                update_u->parent->remove_child(update_u);
                update_u->parent = nullptr;
                root_clusters.push_back(update_u);
            }
        }
        // Delete the level i+1 del clusters that should actually be deleted
        for (auto cluster : del_clusters) {
            add_next_del(cluster->parent);
            if (should_delete(cluster)) {
                disconnect_children(cluster);
                FOR_ALL_NEIGHBORS(cluster, [&](Cluster* neighbor, e_t _) {
                    neighbor->remove_neighbor(cluster);
                });
                if (cluster->parent) cluster->parent->remove_child(cluster);
                free_cluster(cluster);
            } else {
                cluster->mark = false;
            }
        }
        // For insertions add the new edge in level i
        if (update_type == INSERT && update_u && update_v) {
            if constexpr (std::is_same<e_t, empty_t>::value) {
                update_u->insert_neighbor(update_v);
                update_v->insert_neighbor(update_u);
            } else {
                update_u->insert_neighbor_with_value(update_v, link_value);
                update_v->insert_neighbor_with_value(update_u, link_value);
            }
            update_u = update_u->parent;
            update_v = update_v->parent;
        }
        // Recluster the level i root clusters
        recluster_roots();
        // Prepare the next level
        std::swap(root_clusters, next_root_clusters);
        next_root_clusters.clear();
        std::swap(del_clusters, next_del_clusters);
        next_del_clusters.clear();
    }
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::recluster_roots() {
    // Merge deg 3-5 root clusters with all of its deg 1 neighbors
    for (auto cluster : root_clusters) {
        if (!cluster->parent && cluster->degree > 2) [[unlikely]] {
            assert(cluster->degree <= 5);
            auto parent = allocate_cluster();
            if constexpr (!std::is_same<e_t, empty_t>::value) {
                parent->value = identity_v;
            }
            parent->insert_child(cluster);
            cluster->parent = parent;
            next_root_clusters.push_back(parent);
            assert(UFO_NEIGHBOR_MAX >= 3);
            FOR_ALL_NEIGHBORS(cluster, [&](Cluster* neighbor, e_t edge_value) {
                if (neighbor->degree == 1) [[unlikely]] {
                    auto curr = neighbor->parent;
                    while (curr && curr->fanout == 1) {
                        auto next = curr->parent;
                        free_cluster(curr);
                        curr = next;
                    }
                    neighbor->parent = parent;
                    parent->insert_child(neighbor);
                } else if (neighbor->parent) { // Populate new parent's neighbors
                    if constexpr (std::is_same<e_t, empty_t>::value) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    } else {
                        parent->insert_neighbor_with_value(neighbor->parent, edge_value);
                        neighbor->parent->insert_neighbor_with_value(parent, edge_value);
                    }
                }
            });
        }
    }
    // This loop handles all deg 2 and 1 root clusters
    for (auto cluster : root_clusters) {
        // Combine deg 2 root clusters with deg 2 root clusters
        if (!cluster->parent && cluster->degree == 2) [[unlikely]] {
            assert(UFO_NEIGHBOR_MAX >= 2);
            for (int i = 0; i < 2; ++i) {
                auto neighbor = cluster->neighbors[i];
                if (!neighbor->parent && (neighbor->degree == 2)) [[unlikely]] {
                    auto parent = allocate_cluster();
                    cluster->parent = parent;
                    neighbor->parent = parent;
                    parent->insert_child(cluster);
                    parent->insert_child(neighbor);
                    if constexpr (!std::is_same<e_t, empty_t>::value) { // Path query
                        parent->value = f_e(cluster->value, f_e(neighbor->value, cluster->get_edge_value(i)));
                    }
                    next_root_clusters.push_back(parent);
                    for (int i = 0; i < 2; ++i) { // Populate new parent's neighbors
                        if (cluster->neighbors[i]->parent && cluster->neighbors[i]->parent != parent) {
                            if constexpr (std::is_same<e_t, empty_t>::value) {
                                parent->insert_neighbor(cluster->neighbors[i]->parent);
                                cluster->neighbors[i]->parent->insert_neighbor(parent);
                            } else {
                                parent->insert_neighbor_with_value(cluster->neighbors[i]->parent, cluster->get_edge_value(i));
                                cluster->neighbors[i]->parent->insert_neighbor_with_value(parent, cluster->get_edge_value(i));
                            }
                        }
                        if (neighbor->neighbors[i]->parent && neighbor->neighbors[i]->parent != parent) {
                            if constexpr (std::is_same<e_t, empty_t>::value) {
                                parent->insert_neighbor(neighbor->neighbors[i]->parent);
                                neighbor->neighbors[i]->parent->insert_neighbor(parent);
                            } else {
                                parent->insert_neighbor_with_value(neighbor->neighbors[i]->parent, neighbor->get_edge_value(i));
                                neighbor->neighbors[i]->parent->insert_neighbor_with_value(parent, neighbor->get_edge_value(i));
                            }
                        }
                    }
                    break;
                }
            }
            // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
            if (!cluster->parent) [[unlikely]] {
                assert(UFO_NEIGHBOR_MAX >= 2);
                for (int i = 0; i < 2; ++i) {
                    auto neighbor = cluster->neighbors[i];
                    if (neighbor->parent && (neighbor->degree == 1 || neighbor->degree == 2)) [[unlikely]] {
                        if (neighbor->contracts()) continue;
                        cluster->parent = neighbor->parent;
                        neighbor->parent->insert_child(cluster);
                        if constexpr (!std::is_same<e_t, empty_t>::value) { // Path query
                            cluster->parent->value = f_e(cluster->value, f_e(neighbor->value, cluster->get_edge_value(i)));
                        }
                        add_next_del(cluster->parent->parent);
                        auto other_neighbor = cluster->neighbors[!i]; // Popoulate neighbors
                        if (other_neighbor->parent) {
                            if constexpr (std::is_same<e_t, empty_t>::value) {
                                cluster->parent->insert_neighbor(other_neighbor->parent);
                                other_neighbor->parent->insert_neighbor(cluster->parent);
                                // insert_adjacency(cluster->parent, other_neighbor->parent);
                            } else {
                                cluster->parent->insert_neighbor_with_value(other_neighbor->parent, cluster->get_edge_value(!i));
                                other_neighbor->parent->insert_neighbor_with_value(cluster->parent, cluster->get_edge_value(!i));
                                // insert_adjacency(cluster->parent, other_neighbor->parent, cluster->get_edge_value(!i));
                            }
                        }
                        break;
                    }
                }
            }
        // Always combine deg 1 root clusters with its neighboring cluster
        } else if (!cluster->parent && cluster->degree == 1) [[unlikely]] {
            auto neighbor = cluster->neighbors[0];
            if (neighbor->parent) {
                if (neighbor->degree == 2 && neighbor->contracts()) continue;
                cluster->parent = neighbor->parent;
                neighbor->parent->insert_child(cluster);
                add_next_del(cluster->parent->parent);
            } else {
                auto parent = allocate_cluster();
                cluster->parent = parent;
                neighbor->parent = parent;
                parent->insert_child(cluster);
                parent->insert_child(neighbor);
                if constexpr (!std::is_same<e_t, empty_t>::value) { // Path query
                    parent->value = identity_v;
                }
                for (int i = 0; i < 2; ++i) { // Populate new parent's neighbors
                    if (neighbor->neighbors[i] && neighbor->neighbors[i] != cluster && neighbor->neighbors[i]->parent) {
                        if constexpr (std::is_same<e_t, empty_t>::value) {
                            parent->insert_neighbor(neighbor->neighbors[i]->parent);
                            neighbor->neighbors[i]->parent->insert_neighbor(parent);
                        } else {
                            parent->insert_neighbor_with_value(neighbor->neighbors[i]->parent, neighbor->get_edge_value(i));
                            neighbor->neighbors[i]->parent->insert_neighbor_with_value(parent, neighbor->get_edge_value(i));
                        }
                    }
                }
                next_root_clusters.push_back(parent);
            }
        }
    }
    // Add remaining uncombined clusters to the next level
    for (auto cluster : root_clusters) {
        if (!cluster->parent && cluster->degree > 0) [[unlikely]] {
            auto parent = allocate_cluster();
            cluster->parent = parent;
            parent->insert_child(cluster);
            if constexpr (!std::is_same<v_t, empty_t>::value) { // Path query
                parent->value = cluster->value;
            }
            for (int i = 0; i < 2; ++i) { // Populate new parent's neighbors
                if (cluster->neighbors[i] && cluster->neighbors[i]->parent) {
                    if constexpr (std::is_same<e_t, empty_t>::value) {
                        parent->insert_neighbor(cluster->neighbors[i]->parent);
                        cluster->neighbors[i]->parent->insert_neighbor(parent);
                    } else {
                        parent->insert_neighbor_with_value(cluster->neighbors[i]->parent, cluster->get_edge_value(i));
                        cluster->neighbors[i]->parent->insert_neighbor_with_value(parent, cluster->get_edge_value(i));
                    }
                }
            }
            next_root_clusters.push_back(parent);
        }
    }
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
template<typename v_t, typename e_t>
bool UFOTree<v_t, e_t>::connected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}

template<typename v_t, typename e_t>
e_t UFOTree<v_t, e_t>::path_query(vertex_t u, vertex_t v) {
    assert(u < leaves.size() && u >= 0 && v < leaves.size() && v >= 0 && u != v && connected(u, v)); 

    e_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = path_v1 = path_v2 = identity_e;
    Cluster *bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
    bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = nullptr;
    if (leaves[u].degree == 2) {
        bdry_u1 = leaves[u].neighbors[0];
        bdry_u2 = leaves[u].neighbors[1];
    }
    if (leaves[v].degree == 2) {
        bdry_v1 = leaves[v].neighbors[0];
        bdry_v2 = leaves[v].neighbors[1];
    }
    auto curr_u = &leaves[u];
    auto curr_v = &leaves[v];
    while (curr_u->parent != curr_v->parent) {
        // NOTE(ATHARVA): Make this all into one function.
        if (curr_u->degree > 2) {
            if (curr_u->parent->degree == 2) {
                // Superunary to Binary
                bdry_u1 = curr_u->parent->neighbors[0];
                bdry_u2 = curr_u->parent->neighbors[1];
                path_u2 = path_u1;
            }
        } else {
            for (int i = 0; i < 2; i++) {
                auto neighbor = curr_u->neighbors[i];
                if (neighbor && neighbor->parent == curr_u->parent) {
                    if (curr_u->degree == 2) {
                        if (curr_u->parent->degree == 2) {
                            // Binary to Binary
                            if (neighbor == bdry_u1) {
                                path_u1 = f_e(path_u1, f_e(curr_u->get_edge_value(i), neighbor->value));
                                bdry_u2 = bdry_u2->parent;
                                for (int i = 0; i < 2; i++)
                                    if (curr_u->parent->neighbors[i] && curr_u->parent->neighbors[i] != bdry_u2)
                                        bdry_u1 = curr_u->parent->neighbors[i];
                            } else {
                                path_u2 = f_e(path_u2, f_e(curr_u->get_edge_value(i), neighbor->value));
                                bdry_u1 = bdry_u1->parent;
                                for (int i = 0; i < 2; i++)
                                    if (curr_u->parent->neighbors[i] && curr_u->parent->neighbors[i] != bdry_u1)
                                        bdry_u2 = curr_u->parent->neighbors[i];
                            }
                        } else {
                            // Binary to Unary
                            path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
                        }
                    } else {
                        if (curr_u->parent->degree == 2) {
                            // Unary to Binary
                            path_u1 = path_u2 = f_e(path_u1, curr_u->get_edge_value(i));
                            bdry_u1 = curr_u->parent->neighbors[0];
                            bdry_u2 = curr_u->parent->neighbors[1];
                        } else {
                            // Unary to Unary and Unary to Superunary
                            path_u1 = f_e(path_u1, f_e(curr_u->get_edge_value(i), neighbor->value));
                        }
                    }
                    break;
                }
            }
            if (!curr_u->contracts()) {
                if (bdry_u1) bdry_u1 = bdry_u1->parent;
                if (bdry_u2) bdry_u2 = bdry_u2->parent;
            }
        }
        curr_u = curr_u->parent;
        // Same thing for the side of curr_v
        if (curr_v->degree > 2) {
            if (curr_v->parent->degree == 2) {
                // Superunary to Superunary/Binary
                bdry_v1 = curr_v->parent->neighbors[0];
                bdry_v2 = curr_v->parent->neighbors[1];
                path_v2 = path_v1;
            }
        } else {
            for (int i = 0; i < 2; i++) {
                auto neighbor = curr_v->neighbors[i];
                if (neighbor && neighbor->parent == curr_v->parent) {
                    if (curr_v->degree == 2) {
                        if (curr_v->parent->degree == 2) {
                            // Binary to Binary
                            if (neighbor == bdry_v1) {
                                path_v1 = f_e(path_v1, f_e(curr_v->get_edge_value(i), neighbor->value));
                                bdry_v2 = bdry_v2->parent;
                                for (int i = 0; i < 2; i++)
                                    if (curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v2)
                                        bdry_v1 = curr_v->parent->neighbors[i];
                            } else {
                                path_v2 = f_e(path_v2, f_e(curr_v->get_edge_value(i), neighbor->value));
                                bdry_v1 = bdry_v1->parent;
                                for (int i = 0; i < 2; i++)
                                    if (curr_v->parent->neighbors[i] && curr_v->parent->neighbors[i] != bdry_v1)
                                        bdry_v2 = curr_v->parent->neighbors[i];
                            }
                        } else {
                            // Binary to Unary
                            path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                        }
                    } else {
                        if (curr_v->parent->degree == 2) {
                            // Unary to Binary
                            path_v1 = path_v2 = f_e(path_v1, curr_v->get_edge_value(i));
                            bdry_v1 = curr_v->parent->neighbors[0];
                            bdry_v2 = curr_v->parent->neighbors[1];
                        } else {
                            // Unary to Unary and Unary to Superunary
                            path_v1 = f_e(path_v1, f_e(curr_v->get_edge_value(i), neighbor->value));
                        }
                    }
                    break;
                }
            }
            if (!curr_v->contracts()) {
                if (bdry_v1) bdry_v1 = bdry_v1->parent;
                if (bdry_v2) bdry_v2 = bdry_v2->parent;
            }
        }
        curr_v = curr_v->parent;
    }
    // Get the correct path sides when the two vertices meet at the LCA
    e_t total = identity_e;
    if (curr_u->degree == 2)
        total = f_e(total, (curr_v == bdry_u1) ? path_u1 : path_u2);
    else
        total = f_e(total, path_u1);
    if (curr_v->degree == 2)
        total = f_e(total, (curr_u == bdry_v1) ? path_v1 : path_v2);
    else
        total = f_e(total, path_v1);
    // If the LCA contracts them in a star merge, take both edges to the center
    if (curr_u->degree == 1 && curr_v->degree == 1
    && curr_u->neighbors[0] != curr_v) [[unlikely]] {
        total = f_e(total, curr_u->get_edge_value(0));
        total = f_e(total, curr_v->get_edge_value(0));
    }
    // Add the value of the last edge (since they contract one must be deg <= 2)
    else [[likely]] {
        for (int i = 0; i < 2; i++) {
            if (curr_u->neighbors[i] == curr_v) {
                total = f_e(total, curr_u->get_edge_value(i));
                break;
            }
            if (curr_v->neighbors[i] == curr_u) {
                total = f_e(total, curr_v->get_edge_value(i));
                break;
            }
        }
    }
    return total;
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::space() {
    std::unordered_set<Cluster*> visited;
    size_t memory = sizeof(UFOTree<v_t, e_t>);
    for (auto cluster : leaves) {
        memory += cluster.calculate_size();
        auto parent = cluster.parent;
        while (parent != nullptr && visited.count(parent) == 0) {
            memory += parent->calculate_size();
            visited.insert(parent);
            parent = parent->parent;
        }
    }
    return memory;
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::count_nodes() {
    std::unordered_set<Cluster*> visited;
    size_t node_count = 0;
    for(auto cluster : leaves){
        node_count += 1;
        auto parent = cluster.parent;
        while(parent != nullptr && visited.count(parent) == 0){
            node_count += 1;
            visited.insert(parent);
            parent = parent->parent;
        }
    }
    return node_count;
}

template<typename v_t, typename e_t>
size_t UFOTree<v_t, e_t>::get_height() {
    size_t max_height = 0;
    for (vertex_t v = 0; v < leaves.size(); ++v) {
        size_t height = 0;
        Cluster* curr = &leaves[v];
        while (curr) {
            height++;
            curr = curr->parent;
        }
        max_height = std::max(max_height, height);
    }
    return max_height;
}
