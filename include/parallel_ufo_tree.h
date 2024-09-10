#pragma once

#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "sequential_forest.h"


// #define PRINT_DEBUG_INFO

using namespace parlay;

template <typename aug_t>
class ParallelUFOTree {
public:
    using ForestStructureType = SequentialForest;
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k, QueryType q = PATH,
        std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y) -> aug_t {
        return x + y; }, aug_t id = 0, aug_t dval = 0);
    void batch_link(sequence<Edge>& links);
    void batch_cut(sequence<Edge>& cuts);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();

private:
    // Class data and parameters
    sequence<ForestStructureType> forests;
    vertex_t n;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    // Store R and D, root and del clusters
    hashbag<vertex_t> R;
    hashbag<vertex_t> D;
    hashbag<vertex_t> next_R;
    hashbag<vertex_t> next_D;
    // Main batch update functions
    void update_tree(sequence<Edge>& updates, bool deletion);
    void recluster_level(int level, bool deletion, sequence<Edge>& U);
    void set_partners(int level);
    void add_parents(int level);
    void add_adjacency(int level);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k, QueryType q,
std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d)
    : n(n), query_type(q), f(f), identity(id), default_value(d), R(n), D(n), next_R(n), next_D(n) {
    forests.emplace_back();
    sequence<vertex_t> leaves = tabulate(n, [&] (size_t i) { return (vertex_t) i; });
    forests[0].insert_vertices(leaves);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(sequence<Edge>& links) {
    #ifdef PRINT_DEBUG_INFO
    std::cout << std::endl << "[ BATCH LINK ]: ";
    for (Edge e : links) std::cout << e.src << "," << e.dst << " ";
    std::cout << std::endl << std::endl;
    #endif
    update_tree(links, false);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(sequence<Edge>& cuts) {
    #ifdef PRINT_DEBUG_INFO
    std::cout << std::endl << "[ BATCH CUT ]: ";
    for (Edge e : cuts) std::cout << e.src << "," << e.dst << " ";
    std::cout << std::endl << std::endl;
    #endif
    update_tree(cuts, true);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::update_tree(sequence<Edge>& updates, bool deletion) {
    // Determine the intial set of root clusters and deletion clusters
    parallel_for(0, updates.size(), [&](size_t i) {
        if (forests[0].try_set_status_atomic(updates[i].src, ROOT)) R.insert(updates[i].src);
        if (forests[0].try_set_status_atomic(updates[i].dst, ROOT)) R.insert(updates[i].dst);
        if (forests.size() > 1) {
            vertex_t parent1 = forests[0].get_parent(updates[i].src);
            vertex_t parent2 = forests[0].get_parent(updates[i].dst);
            if (parent1 != NONE && forests[1].try_set_status_atomic(parent1, DEL)) D.insert(parent1);
            if (parent2 != NONE && forests[1].try_set_status_atomic(parent2, DEL)) D.insert(parent2);
        }
    });
    // Start the level-by-level update procedure
    recluster_level(0, deletion, updates);
    #ifdef PRINT_DEBUG_INFO
    print_tree();
    #endif
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_level(int level, bool deletion, sequence<Edge>& U) {
    #ifdef PRINT_DEBUG_INFO
    std::cout << "RECLUSTER LEVEL " << level << std::endl;
    std::cout << "R ";
    R.print();
    std::cout << "D ";
    D.print();
    #endif
    bool next_R_empty = true;
    bool next_D_empty = true;
    // Remove the parents of root clusters whose parent will get deleted/unset
    sequence<vertex_t> disconnect_from_parent = parlay::map_maybe(R.extract_all(), [&] (auto v) -> std::optional<vertex_t> {
        vertex_t parent = forests[level].get_parent(v);
        if (parent != NONE) {
            bool low_degree = forests[level+1].get_degree(parent) < 3;
            bool low_fanout = forests[level+1].get_child_count(parent) < 3;
            if (low_degree && low_fanout || forests[level].get_degree(v) < 3) {
                forests[level].mark(v);
                return forests[level].get_parent(v);
            }
        }
        return std::nullopt;
    });
    if (forests.size() > level+1) forests[level+1].subtract_children(disconnect_from_parent);
    // At level 0 remove root clusters whose parent is not unset
    if (level == 0) {
        R.for_all([&](vertex_t v) {
            if (forests[level].get_parent(v) == NONE || forests[level].is_marked(v)) next_R.insert(v);
            else forests[level].unset_status(v);
        });
        R.clear();
        std::swap(R, next_R);
    }
    // Delete the edges from our initial updates that are still in level i+1
    if (deletion) {
        auto next_U = forests[level].map_edges_to_parents(U);
        if (forests.size() > level+1) forests[level+1].delete_edges(next_U);
    }

    // Collect the new level i root clusters that will be formed by deleting some level i+1 clusters
    // Also mark root clusters whose parents should be deleted/unset
    R.for_all([&](vertex_t v) {
        vertex_t parent = forests[level].get_parent(v);
        if (parent != NONE) {
            bool low_degree = forests[level+1].get_degree(parent) < 3;
            bool low_fanout = forests[level+1].get_child_count(parent) < 3;
            if (low_degree && low_fanout) {
                auto iter = forests[level].get_neighbor_iterator(v);
                for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                    if (forests[level].get_parent(neighbor) == parent) {
                        if (forests[level].try_set_status_atomic(neighbor, ROOT)) {
                            R.insert(neighbor);
                            forests[level].mark(neighbor);
                            auto iter = forests[level].get_neighbor_iterator(neighbor);
                            for(vertex_t grand_neighbor = iter->next(); grand_neighbor != NONE; grand_neighbor = iter->next()) {
                                if (forests[level].get_parent(grand_neighbor) == parent) {
                                    if (forests[level].try_set_status_atomic(grand_neighbor, ROOT)) {
                                        R.insert(grand_neighbor);
                                        forests[level].mark(grand_neighbor);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    });

    // Delete the level 0 edges in the initial batch
    if (level == 0) {
        if (deletion) {
            forests[0].delete_edges(U);
        } else {
            forests[0].insert_edges(U);
        }
        U = forests[0].map_edges_to_parents(U);
    }
    if (forests.size() > level+1) U = forests[level+1].map_edges_to_parents(U);

    // Unset the parents of all marked clusters
    R.for_all([&](vertex_t v) {
        if (forests[level].is_marked(v)) {
            forests[level].unset_parent(v);
            forests[level].unmark(v);
        }
    });

    // Determine the new combinations of the root clusters at this level
    set_partners(level);

    // Get the next set of del clusters and do the deletion of clusters
    sequence<vertex_t> del = filter(D.extract_all(), [&] (auto v) {
        vertex_t parent = forests[level+1].get_parent(v);
        bool low_degree = forests[level+1].get_degree(v) < 3;
        bool low_fanout = forests[level+1].get_child_count(v) < 3;
        if (!low_degree || !low_fanout) { // Current del clusters that were not deleted added to next_R
            forests[level+1].set_status(v, ROOT);
            next_R.insert(v);
            next_R_empty = false;
        }
        if (parent != NONE) {
            if (forests[level+2].try_set_status_atomic(parent, DEL)) { // Parents added to next_D
                next_D.insert(parent);
                next_D_empty = false;
            }
        }
        return (low_degree && low_fanout);
    });
    sequence<vertex_t> del_parents = parlay::map(del, [&] (auto v) { return forests[level+1].get_parent(v); });
    if (forests.size() > level+2) forests[level+2].subtract_children(del_parents);
    D.for_all([&](vertex_t v) { // Non-deleted neighbors of parents that will be deleted added to next_R
        vertex_t parent = forests[level+1].get_parent(v);
        if (parent != NONE) {
            bool p_low_degree = forests[level+2].get_degree(parent) < 3;
            bool p_low_fanout = forests[level+2].get_child_count(parent) < 3;
            if (p_low_degree && p_low_fanout) {
                auto iter = forests[level+1].get_neighbor_iterator(v);
                for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                    if (forests[level+1].get_parent(neighbor) == parent) {
                        if (forests[level+1].try_set_status_atomic(neighbor, ROOT)) {
                            next_R.insert(neighbor);
                            next_R_empty = false;
                            auto iter = forests[level+1].get_neighbor_iterator(neighbor);
                            for(vertex_t grand_neighbor = iter->next(); grand_neighbor != NONE; grand_neighbor = iter->next()) {
                                if (forests[level+1].get_parent(grand_neighbor) == parent) {
                                    if (forests[level+1].try_set_status_atomic(grand_neighbor, ROOT)) {
                                        next_R.insert(grand_neighbor);
                                        next_R_empty = false;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    });
    if (forests.size() > level+1) forests[level+1].delete_vertices(del);

    // Add the new parents and edges to the next level for the set of combinations determined
    add_parents(level);
    add_adjacency(level);

    // Prepare the inputs to the next level
    R.for_all([&](vertex_t v) {
        vertex_t parent = forests[level].get_parent(v);
        if (parent != NONE) {
            if (forests[level+1].try_set_status_atomic(parent, ROOT)) {
                next_R.insert(parent);
                next_R_empty = false;
            }
            vertex_t p_parent = forests[level+1].get_parent(parent);
            if (p_parent != NONE) {
                if (forests[level+2].try_set_status_atomic(p_parent, DEL)) {
                    next_D.insert(p_parent);
                    next_D_empty = false;
                }
            }
        }
    });
    // Cleanup
    R.for_all([&](vertex_t v) {
        vertex_t partner = forests[level].get_partner(v);
        if (partner != NONE) forests[level].unset_partner(partner);
        forests[level].unset_partner(v);
        forests[level].unset_status(v);
    });
    R.clear();
    D.clear();
    // Check if the update is complete or recursively recluster the next level
    if (next_R_empty && next_D_empty) return;
    std::swap(R, next_R);
    std::swap(D, next_D);
    recluster_level(level+1, deletion, U);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::set_partners(int level) {
    // Perform clustering of the root clusters
    R.for_all([&](vertex_t cluster) {
        if (forests[level].get_parent(cluster) != NONE) return;
        // Combine deg 3+ root clusters with deg 1 root or non-root clusters
        if (forests[level].get_degree(cluster) >= 3) {
            // If it became high degree, find any non-contracting deg 1 neighbors
            auto iter = forests[level].get_neighbor_iterator(cluster);
            for (vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                if (forests[level].get_degree(neighbor) == 1) {
                    forests[level].set_partner(cluster, neighbor);
                    forests[level].set_partner(neighbor, cluster);
                    vertex_t parent = forests[level].get_parent(neighbor);
                    if (parent != NONE) {
                        if (forests[level+1].try_set_status_atomic(parent, DEL)) {
                            forests[level].unset_parent(neighbor);
                            R.insert(neighbor);
                            D.insert(parent);
                        }
                    }
                }
            }
        } else if (forests[level].get_degree(cluster) == 2) {
            // Only local maxima in priority with respect to deg 2 clusters will act
            if (!forests[level].is_local_max_priority(cluster)) return;
            // Travel left/right and pair clusters until a deg 3, deg 1, non-root,
            // or combined cluster found
            vertex_t first_neighbor = forests[level].get_first_neighbor(cluster);
            for (auto direction : {0, 1}) {
                auto curr = cluster;
                auto next = direction ? first_neighbor : forests[level].get_other_neighbor(cluster, first_neighbor);
                if (forests[level].get_partner(curr) != NONE) {
                    curr = next;
                    next = forests[level].get_other_neighbor(curr, cluster);
                }
                while (curr != NONE && forests[level].get_parent(curr) == NONE && forests[level].get_degree(curr) == 2
                && next != NONE && forests[level].get_degree(next) < 3 && !forests[level].contracts(next)) {
                    if (!forests[level].try_set_partner_atomic(curr, next)) break;
                    if (forests[level].get_degree(next) == 1) { // If next deg 1 they can combine
                        forests[level].set_partner(next, curr);
                        break;
                    }
                    if (!forests[level].try_set_partner_atomic(next, curr)) { // If the CAS fails next was combined from the other side
                        if (forests[level].get_partner(next) != curr) forests[level].unset_partner(curr);
                        break;
                    }
                    if (forests[level].get_parent(next) != NONE) break; // Stop traversing at a non-root cluster
                    // Get the next two clusters in the chain
                    curr = forests[level].get_other_neighbor(next, curr);
                    if (curr != NONE) next = forests[level].get_other_neighbor(curr, next);
                }
            }
        } else if (forests[level].get_degree(cluster) == 1) {
            auto iter = forests[level].get_neighbor_iterator(cluster);
            for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                if (forests[level].get_degree(neighbor) == 1) {
                    forests[level].set_partner(cluster, neighbor);
                    forests[level].set_partner(neighbor, cluster);
                    break;
                }
                // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
                if (forests[level].get_parent(neighbor) != NONE && forests[level].get_degree(neighbor) == 2) {
                    if (forests[level].contracts(neighbor)) continue;
                    if (!forests[level].try_set_partner_atomic(neighbor, cluster)) continue;
                    forests[level].set_partner(cluster, neighbor);
                    break;
                }
                // Combine deg 1 root cluster with deg 3+ clusters always
                if (forests[level].get_degree(neighbor) >= 3) {
                    forests[level].set_partner(cluster, neighbor);
                    forests[level].set_partner(neighbor, cluster);
                    // QUINTEN: if a high degree cluster is a non-root cluster, it must have been high degree previously ?
                    // If it became high degree, any non-contracting deg 1 clusters must be in the first 2 neighbors
                    // if (forests[level].get_parent(neighbor) != NONE) {
                    //     auto n_iter = forests[level].get_neighbor_iterator(neighbor);
                    //     while (vertex_t entry = n_iter->next() != NONE) {
                    //         if (forests[level].get_degree(entry) == 1 && forests[level].get_parent(entry) != forests[level].get_parent(neighbor)) {
                    //             forests[level].try_set_partner(entry, neighbor);
                    //             forests[level].set_parent(entry, forests[level].get_parent(neighbor));
                    //         }
                    //     }
                    // }
                    break;
                }
            }
        }
    });
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::add_parents(int level) {
    // Create new parent clusters for each new combination
    R.for_all([&](vertex_t cluster) {
        vertex_t partner = forests[level].get_partner(cluster);
        vertex_t parent = forests[level].get_parent(cluster);
        // Deal with possible multi-combinations
        if (partner != NONE) {
            vertex_t p_partner = forests[level].get_partner(partner);
            vertex_t p_parent = forests[level].get_parent(partner);
            if (forests[level].get_degree(partner) >= 3) {
                if (p_parent != NONE) forests[level].set_parent(cluster, p_parent);
                return;
            }
            if (forests[level].get_degree(cluster) >= 3 && parent == NONE) {
                auto new_parent = cluster;
                forests[level].set_parent(cluster, new_parent);
                auto iter = forests[level].get_neighbor_iterator(cluster);
                for (vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                    if (forests[level].get_partner(neighbor) == cluster) {
                        forests[level].set_parent(neighbor, new_parent);
                    }
                }
                return;
            }
        }
        // Deal with regular pairwise combinations
        if (partner != NONE) {
            vertex_t p_partner = forests[level].get_partner(partner);
            vertex_t p_parent = forests[level].get_parent(partner);
            if (p_parent != NONE && p_parent != parent) {   // partner is non-root cluster
                forests[level].set_parent(cluster, p_parent);
                forests[level].unset_partner(partner);
        } else {                                            // two root clusters
                vertex_t new_parent = std::max(cluster, partner);
                forests[level].set_parent(cluster, new_parent);
                forests[level].set_parent(partner, new_parent);
            }
        } else if (parent == NONE && forests[level].get_degree(cluster) > 0) { // clusters that don't combine get a new parent
            vertex_t new_parent = cluster;
            forests[level].set_parent(cluster, new_parent);
            forests[level].set_partner(cluster, cluster);
        }
    });
    // Add the parents
    sequence<vertex_t> P = parlay::map(R.extract_all(), [&](vertex_t v) {
        return forests[level].get_parent(v);
    });
    if (forests.size() <= level+1) forests.emplace_back();
    forests[level+1].insert_vertices(P);
    forests[level+1].add_children(P);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::add_adjacency(int level) {
    // Fill in the neighbor lists of the new clusters
    hashbag<edge_t> E(n);
    R.for_all([&](vertex_t v) {
        auto iter = forests[level].get_neighbor_iterator(v);
        for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
            E.insert(VERTICES_TO_EDGE(v, neighbor));
        }
    });
    sequence<Edge> edges = parlay::map_maybe(E.extract_all(), [&](edge_t e) -> std::optional<Edge> {
        Edge edge = EDGE_TYPE_TO_STRUCT(e);
        if (forests[level].get_parent(edge.src) != forests[level].get_parent(edge.dst))
            return (Edge) {forests[level].get_parent(edge.src), forests[level].get_parent(edge.dst)};
        return std::nullopt;
    });
    if (forests.size() > level+1) forests[level+1].insert_edges(edges);
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    vertex_t root_u = u;
    int level = 0;
    while (forests[level].get_parent(root_u) != -1)
        root_u = forests[level++].get_parent(root_u);
    vertex_t root_v = v;
    level = 0;
    while (forests[level].get_parent(root_v) != -1)
        root_v = forests[level++].get_parent(root_v);
    return root_u == root_v;
}
