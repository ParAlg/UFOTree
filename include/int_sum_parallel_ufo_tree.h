#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include "aggregator.h"
#include "int_sum_parallel_ufo_cluster.h"


namespace ufo::int_sum {

template <typename aug_t = empty_t>
class ParallelUFOTree {
    using Cluster = ParallelUFOCluster<aug_t>;
    static constexpr uintptr_t NULL_PTR = 0;
    static constexpr uintptr_t NULL_PAR = 1;
    static constexpr uintptr_t DEL_MARK = 2;
    static constexpr uintptr_t NON_ROOT_MARK = 3;
    static constexpr uintptr_t NEW_PAR_MARK = 4;
    using allocator = parlay::type_allocator<Cluster>;
    using EdgeSlice = parlay::slice<std::pair<Cluster*, Cluster*>*, std::pair<Cluster*, Cluster*>*>;
    static inline uint32_t combine_weights(uint32_t a, uint32_t b) { return a + b; }
    static inline uint32_t recompute_parent_aggregate_from_pair(Cluster* c1, Cluster* c2) {
        if (!c1 || !c2) return 0;
        int d1 = c1->get_degree();
        int d2 = c2->get_degree();
        if (d1 == 2 && d2 == 2) return combine_weights(c1->aggregate_weight, c2->aggregate_weight);
        if (d1 >= 3 && d2 == 1) return c1->aggregate_weight;
        if (d2 >= 3 && d1 == 1) return c2->aggregate_weight;
        return 0;
    }
    static inline uint32_t recompute_parent_aggregate_from_child(Cluster* child) {
        Cluster* parent = child->parent;
        if (!parent) return child->aggregate_weight;

        int child_degree = child->get_degree();
        if (child_degree >= 3) return child->aggregate_weight;

        if (!child->contracts()) return child->aggregate_weight;

        if (child_degree == 2) {
            Cluster* n1 = child->get_neighbor();
            Cluster* n2 = child->get_other_neighbor(n1);
            Cluster* partner = nullptr;
            if (n1 && n1->parent == parent) partner = n1;
            else if (n2 && n2->parent == parent) partner = n2;
            if (!partner) return child->aggregate_weight;
            return recompute_parent_aggregate_from_pair(child, partner);
        }

        if (child_degree == 1) {
            Cluster* partner = child->get_neighbor();
            if (partner && partner->parent == parent) return recompute_parent_aggregate_from_pair(child, partner);
            return child->aggregate_weight;
        }

        return child->aggregate_weight;
    }
    static inline Cluster* advance_to_nonzero_cluster(Cluster* c) {
        // Some transitional internal nodes can carry zero aggregate; walk through unary chains.
        for (int i = 0; i < 8 && c && c->aggregate_weight == 0 && c->get_degree() == 1; ++i) {
            Cluster* next = c->get_neighbor();
            if (!next || next == c) break;
            c = next;
        }
        return c;
    }
    static inline bool is_endpoint_leaf_parent_duplicate(const std::vector<Cluster>& leaf_clusters,
                                                         Cluster* c, vertex_t u, vertex_t v) {
        if (!c) return false;
        for (vertex_t x : {u, v}) {
            if (leaf_clusters[x].parent == c && c->aggregate_weight == leaf_clusters[x].aggregate_weight)
                return true;
        }
        return false;
    }
    static inline bool is_endpoint_subtree_aggregate(const std::vector<Cluster>& leaf_clusters,
                                                     Cluster* c, vertex_t u, vertex_t v) {
        if (!c || c->aggregate_weight == 0) return false;
        for (vertex_t x : {u, v}) {
            if (c->aggregate_weight != leaf_clusters[x].aggregate_weight) continue;
            for (Cluster* p = leaf_clusters[x].parent; p; p = p->parent) {
                if (p == c) return true;
            }
        }
        return false;
    }
public:
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k);
    ~ParallelUFOTree();
    void batch_link(parlay::sequence<std::pair<int, int>>& links);
    void batch_cut(parlay::sequence<std::pair<int, int>>& cuts);
    void BatchLink(parlay::sequence<std::pair<int, int>>& links) { batch_link(links); };
    void BatchCut(parlay::sequence<std::pair<int, int>>& cuts) { batch_cut(cuts); };
    void UpdateWeight(vertex_t v, uint32_t weight);
    uint32_t PathQuery(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid(parlay::sequence<std::pair<int, int>>& edges);
    int get_height(vertex_t v);
    void print_tree();
// private:
    // Class data and parameters
    std::vector<Cluster> leaves;

    // Thread local root cluster aggregator
    Aggregator<Cluster*>* thread_local_root_clusters;
    Aggregator<Cluster*>* thread_local_next_root_clusters;
    Aggregator<std::pair<Cluster*, Cluster*>> thread_local_del_clusters;
    Aggregator<std::pair<Cluster*, Cluster*>> thread_local_dir_edges;
    Aggregator<std::pair<Cluster*, Cluster*>> thread_local_del_edges;
    Aggregator<std::pair<Cluster*, Cluster*>> thread_local_new_del_edges;

    // Helper functions
    void recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type);
    parlay::sequence<std::pair<Cluster*, Cluster*>> process_initial_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups);
    parlay::sequence<std::pair<Cluster*, Cluster*>> process_del_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups);
    void recluster_root_clusters(UpdateType update_type);
    inline void recluster_degree_one_root(Cluster* root_cluster, UpdateType update_type);
    inline void recluster_degree_two_root(Cluster* root_cluster);
    inline void recluster_high_degree_root(Cluster* root_cluster);
    static inline bool is_local_max(Cluster* c);
    void create_new_parents();
    void finish_reclustering();

    static Cluster* allocate_cluster();
    static void free_cluster(Cluster* c);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : leaves(n), thread_local_root_clusters(new Aggregator<Cluster*>()),
    thread_local_next_root_clusters(new Aggregator<Cluster*>()), thread_local_del_clusters(), thread_local_dir_edges(),
    thread_local_del_edges(), thread_local_new_del_edges() {
    // Start with unit vertex weights by default.
    parlay::parallel_for(0, n, [&] (size_t i) {
        leaves[i].aggregate_weight = 1;
    });
}

namespace ufo {
template <typename aug_t = empty_t>
using IntSumParallelUFOTree = int_sum::ParallelUFOTree<aug_t>;
}

template <typename aug_t>
ParallelUFOTree<aug_t>::~ParallelUFOTree() {
    auto clusters_to_delete = parlay::flatten(parlay::tabulate(leaves.size(), [&] (size_t i) {
        parlay::sequence<Cluster*> clusters;
        Cluster* curr = leaves[i].parent;
        while (curr && curr != (Cluster*) DEL_MARK) {
            Cluster* next = AtomicLoad(&curr->parent);
            if (next != (Cluster*) DEL_MARK && CAS(&curr->parent, next, (Cluster*) DEL_MARK))
                clusters.push_back(curr);
            else break;
            curr = next;
        }
        return clusters;
    }));
    parlay::parallel_for(0, clusters_to_delete.size(), [&] (size_t i) {
        free_cluster(clusters_to_delete[i]);
    });
    allocator::finish();
    delete thread_local_root_clusters;
    delete thread_local_next_root_clusters;
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(parlay::sequence<std::pair<int, int>>& links) {
    recluster_tree(links, INSERT);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(parlay::sequence<std::pair<int, int>>& cuts) {
    recluster_tree(cuts, DELETE);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type) {
    // ==============
    // INITIALIZATION
    // ==============
    parlay::sequence<std::pair<Cluster*, Cluster*>> del_clusters;
    parlay::sequence<std::pair<Cluster*, Cluster*>> next_del_clusters;
    parlay::sequence<parlay::sequence<std::pair<Cluster*, Cluster*>>> all_del_clusters;

    // The intial dir updates are just the batch of links or cuts in both directions.
    // For deletion batches, the initial level 1 del edges are initial cuts one level up.
    auto dir_updates = parlay::tabulate(2*updates.size(), [&] (size_t i) {
        Cluster* c1 = &leaves[updates[i/2].first];
        Cluster* c2 = &leaves[updates[i/2].second];
        if (update_type == DELETE) {
            if (i % 2 == 0 && c1->parent && c2->parent && c1->parent != c2->parent) {
                thread_local_del_edges.push_back(std::make_pair(c1->parent, c2->parent));
                thread_local_del_edges.push_back(std::make_pair(c2->parent, c1->parent));
            }
        }
        if (i % 2 == 0) return std::make_pair(c1, c2);
        return std::make_pair(c2, c1);
    });
    auto dir_update_groups = group_by_key_inplace(dir_updates);

    parlay::sequence<std::pair<Cluster*, Cluster*>> parents;
    parlay::sequence<std::pair<Cluster*, EdgeSlice>> parent_groups;
    parlay::parallel_do(
    [&] () {
        // Group the affected vertices by parent.
        parents = parlay::map(dir_update_groups, [&] (auto group) {
            return std::make_pair(group.first->parent, group.first);
        });
        parent_groups = group_by_key_inplace(parents);
    },
    [&] () {
        // For deletion batches, keep a trace of the level i+2 del edges to delete and decrement `degree` at each level.
        if (update_type == DELETE) {
            auto del_edges = thread_local_del_edges.to_sequence();
            thread_local_del_edges.clear();
            auto del_edge_groups = group_by_key_inplace(del_edges);
            parlay::parallel_for(0, del_edge_groups.size(), [&] (size_t i) {
                auto& [cluster, edges] = del_edge_groups[i];
                cluster->degree = cluster->get_degree() - edges.size();
                parlay::parallel_for(0, edges.size(), [&] (size_t j) {
                    if (edges[j].first->parent && edges[j].second->parent && edges[j].first->parent != edges[j].second->parent)
                        thread_local_del_edges.push_back(std::make_pair(edges[j].first->parent, edges[j].second->parent));
                });
            });
        }
    });

    // The initial root clusters are children of level 1 parents that will get deleted, or deg 1 endpoints.
    // The intial del clusters are just all parents of vertices that got an update.
    // This function adds level 0 root clusters to `thread_local_root_clusters`, and returns the level 1 del clusters.
    del_clusters = process_initial_clusters(parent_groups);


    /* This is the main loop over all levels of the tree. Before each level i, `root_clusters` contains the
    level i clusters formed from contractions in the previous level, along with any additional children of level
    i+1 del clusters that will get deleted. All root clusters have their parent pointer cleared already.
    The sequence `del_clusters` contains the level i+1 clusters that were parents of del clusters at the previous
    level, along with the parents of any nodes that started a new remove ancestor path due to a root cluster
    merging into a non-root cluster at the previous level. The sequence 'dir_updates' contains the list of the
    intial update edges that still exist in level i.
    (PHASE 1) First we update level i with 'dir_updates'. At level 0 this corresponds to updating the level 0
    forest. Since some clusters never get deleted, it is necessary to propagate this change to all levels.
    This must happen before the reclustering at level i, but after the del clusters are processed in phase 3,
    so that the recluster has accurate tree information, but phase 3 retains the structure of the previous tree.
    (PHASE 2) Next we, will determine the new set of contractions over the level i root clusters using an MIS. We
    then allocate new level i+1 parent nodes and update the parent pointers of the root clusters.
    These new parents become the `next_root_clusters` along with those determined phase 3. Any contraction
    with a non-root cluster adds the level i+1 parent to `del_clusters`.
    (PHASE 3) First we take the level i+1 `del_clusters`, and populate `next_del_clusters` by mapping them to their
    level i+2 parents. We also determine which level i+2 del clusters will actually be deleted, and mark those, not
    yet deleting them. For deletion batches we first need to update the `degree` field in each level i+2 del cluster
    before deciding if it should be deleted. To do this, we keep track of the original cut edges still alive at this
    level, as well as edges that get deleted incident to lower level deleted clusters. During this phase we delete
    the level i+1 del clusters that were marked at the previous level. We also find any level i+1 non-del children
    of the level i+2 clusters which will be deleted. These will become `next_root_clusters` at the next level i+1.
    (PHASE 4) In this phase we insert any edge incident to a level i root cluster into level i+1, if it wasn't
    contracted away. This fills the adjacency lists of our new parent clusters and also updates any level i+1
    clusters that remain from the previous tree. This must happen after phase 3 so that we do not modify the
    structure of the old tree before processing the del clusters. */
    while (!thread_local_root_clusters->empty() || (del_clusters.size() > 0 && del_clusters[0].first != (Cluster*) NULL_PAR)) {
        // =======
        // PHASE 1
        // =======

        // For insertion batches, insert the linked edges at level i, and map to the next level.
        parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = dir_update_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            if (update_type == INSERT) cluster->insert_neighbors_sorted(neighbors);
            else cluster->delete_neighbors_sorted(neighbors);
            parlay::parallel_for(0, edges.size(), [&] (size_t j) {
                if (edges[j].first->parent && edges[j].second->parent && edges[j].first->parent != edges[j].second->parent)
                    thread_local_dir_edges.push_back(std::make_pair(edges[j].first->parent, edges[j].second->parent));
            });
        });
        dir_updates = thread_local_dir_edges.to_sequence();
        thread_local_dir_edges.clear();
        dir_update_groups = group_by_key_inplace(dir_updates);

        // =======
        // PHASE 2
        // =======

        // Recluster the root clusters.
        recluster_root_clusters(update_type);
        del_clusters = parlay::append(del_clusters, thread_local_del_clusters.to_sequence());
        thread_local_del_clusters.clear();

        parlay::parallel_do(
        [&] () {
            // This returns only the new clusters that were created during the reclustering at this level.
            create_new_parents();
        },
        [&] () {
            // =======
            // PHASE 3
            // =======

            // Determine pointers to level i+1 del clusters that will get deleted.
            // For deletion batches, map deleting edges to level i+2 and add them to the del edges.
            parlay::parallel_for(0, del_clusters.size(), [&] (size_t i) {
                auto [parent, cluster] = del_clusters[i];
                if (parent != (Cluster*) NULL_PTR && parent != (Cluster*) NULL_PAR && cluster->partner == (Cluster*) DEL_MARK) {
                    cluster->for_all_neighbors([&] (Cluster* neighbor) {
                        if (!neighbor->partner) {
                            thread_local_new_del_edges.push_back(std::make_pair(neighbor, cluster));
                            if (update_type == DELETE)
                                if (neighbor->parent && cluster->parent && neighbor->parent != cluster->parent)
                                    thread_local_del_edges.push_back(std::make_pair(neighbor->parent, cluster->parent));
                        }
                    });
                }
            });
            auto new_del_edges = thread_local_new_del_edges.to_sequence();
            thread_local_new_del_edges.clear();
            parlay::sequence<std::pair<Cluster*, EdgeSlice>> new_del_edge_groups;

            parlay::parallel_do(
            [&] () {
                // Group edges to delete by endpoint.
                new_del_edge_groups = group_by_key_inplace(new_del_edges);
            },
            [&] () {
                parlay::parallel_do(
                [&] () {
                    // Group del clusters by parent.
                    parent_groups = group_by_key_inplace(del_clusters);
                },
                [&] () {
                    // For deletion batches, keep a trace of the level i+2 del edges to delete and decrement `degree` at each level.
                    if (update_type == DELETE) {
                        // Group the level i+2 updates by endpoint.
                        auto del_edges = thread_local_del_edges.to_sequence();
                        thread_local_del_edges.clear();
                        auto del_edge_groups = group_by_key_inplace(del_edges);
                        // Update the degree fields and map the level i+2 del edges to level i+3.
                        parlay::parallel_for(0, del_edge_groups.size(), [&] (size_t i) {
                            auto& [cluster, edges] = del_edge_groups[i];
                            cluster->degree = cluster->get_degree() - parlay::unique(edges).size();
                            parlay::parallel_for(0, edges.size(), [&] (size_t j) {
                                if (edges[j].first->parent && edges[j].second->parent && edges[j].first->parent != edges[j].second->parent)
                                    thread_local_del_edges.push_back(std::make_pair(edges[j].first->parent, edges[j].second->parent));
                            });
                        });
                    }
                });
            });

            // The next level i+1 root clusters, which are children of a level i+2 del cluster that will be deleted.
            // This function populates `thread_local_next_root_clusters` with the next level i+1 root clusters.
            // This function returns the next level i+2 del clusters and we store it in `next_del_clusters`.
            next_del_clusters = process_del_clusters(parent_groups);

            // Delete pointers to level i+1 del clusters that will get deleted.
            parlay::parallel_for(0, new_del_edge_groups.size(), [&] (size_t i) {
                auto& [cluster, edges] = new_del_edge_groups[i];
                auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
                cluster->delete_neighbors_sorted(neighbors);
            });

            // Delete level i+1 del clusters that should be deleted.
            all_del_clusters.push_back(std::move(del_clusters));
        });

        // =======
        // PHASE 4
        // =======

        // Clear the partner fields and populate neighbors of level i+1 clusters.
        finish_reclustering();

        // ==================
        // PREPARE NEXT LEVEL
        // ==================
        std::swap(thread_local_root_clusters, thread_local_next_root_clusters);
        thread_local_next_root_clusters->clear();
        del_clusters = std::move(next_del_clusters);
    }

    // Delete all the del_clusters throughout the entire update.
    parlay::parallel_for(0, all_del_clusters.size(), [&] (size_t i) {
        parlay::parallel_for(0, all_del_clusters[i].size(), [&] (size_t j) {
            if (all_del_clusters[i][j].second)
                if (all_del_clusters[i][j].second->partner)
                    free_cluster(all_del_clusters[i][j].second);
        });
    });
}

// ================================================================================================
// ================================        HELPER FUNCTIONS        ================================
// ================================================================================================

template <typename aug_t>
parlay::sequence<std::pair<ParallelUFOCluster<aug_t>*, ParallelUFOCluster<aug_t>*>> ParallelUFOTree<aug_t>::process_initial_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups) {
    return parlay::tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        if (parent == (Cluster*) NULL_PTR) {
            parlay::parallel_for(0, children.size(), [&] (size_t i) {
                thread_local_root_clusters->push_back(children[i].second);
            });
            return std::make_pair((Cluster*) NULL_PAR, (Cluster*) NULL_PTR);
        }

        Cluster* max;
        int max_degree;
        if (children.size() > 500) {
            max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
            max_degree = max->get_degree();
        } else {
            max = children[0].second;
            max_degree = max->get_degree();
            for (size_t j = 1; j < children.size(); ++j) {
                if (children[j].second->get_degree() > max_degree) {
                    max = children[j].second;
                    max_degree = max->get_degree();
                }
            }
        }
        Cluster* center = max;

        if (max_degree == 1) {
            center = max->get_neighbor();
            if (AtomicLoad(&center->parent) != parent) {
                thread_local_root_clusters->push_back(max);
                AtomicStore(&max->parent, (Cluster*) NULL_PTR);
                parent->partner = (Cluster*) DEL_MARK;
                return std::make_pair(parent->parent, parent);
            }
        }

        int fanout = center->get_degree() - parent->get_degree() - children.size();
        if (center == max) fanout++;
        int degree = parent->degree;

        if (fanout < 4 && degree < 4) {
            center->for_all_neighbors([&] (auto neighbor) {
                if(AtomicLoad(&neighbor->parent) == parent) {
                    thread_local_root_clusters->push_back(neighbor);
                    AtomicStore(&neighbor->parent, (Cluster*) NULL_PTR);
                }
            });
            thread_local_root_clusters->push_back(center);
            AtomicStore(&center->parent, (Cluster*) NULL_PTR);
            parent->partner = (Cluster*) DEL_MARK;
        } else {
            if (center == max) {
                if (children.size() > 500) {
                    parlay::parallel_for(0, children.size(), [&] (size_t i) {
                        if (children[i].second != max) {
                            thread_local_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    });
                } else {
                    for (size_t i = 0; i < children.size(); ++i) {
                        if (children[i].second != max) {
                            thread_local_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    }
                }
            } else {
                if (children.size() > 500) {
                    parlay::parallel_for(0, children.size(), [&] (size_t i) {
                        thread_local_root_clusters->push_back(children[i].second);
                        AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                    });
                } else {
                    for (size_t i = 0; i < children.size(); ++i) {
                        thread_local_root_clusters->push_back(children[i].second);
                        AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                    }
                }
            }
        }

        return std::make_pair(parent->parent, parent);
    });
}

template <typename aug_t>
parlay::sequence<std::pair<ParallelUFOCluster<aug_t>*, ParallelUFOCluster<aug_t>*>> ParallelUFOTree<aug_t>::process_del_clusters(parlay::sequence<std::pair<Cluster*, EdgeSlice>>& parent_groups) {
    return parlay::tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        if (parent == (Cluster*) NULL_PAR) return std::make_pair((Cluster*) NULL_PAR, (Cluster*) NULL_PTR);
        if (parent == (Cluster*) NULL_PTR) {
            parlay::parallel_for(0, children.size(), [&] (size_t i) {
                if (!children[i].second->partner) thread_local_next_root_clusters->push_back(children[i].second);
            });
            return std::make_pair((Cluster*) NULL_PAR, (Cluster*) NULL_PTR);
        }

        Cluster* max;
        int max_degree;
        if (children.size() > 500) {
            max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
            max_degree = max->get_degree();
        } else {
            max = children[0].second;
            max_degree = max->get_degree();
            for (size_t j = 1; j < children.size(); ++j) {
                if (children[j].second->get_degree() > max_degree) {
                    max = children[j].second;
                    max_degree = max->get_degree();
                }
            }
        }
        Cluster* center = max;

        if (max_degree == 1) {
            center = max->get_neighbor();
            if (AtomicLoad(&center->parent) != parent) {
                if (!max->partner) {
                    thread_local_next_root_clusters->push_back(max);
                    AtomicStore(&max->parent, (Cluster*) NULL_PTR);
                }
                parent->partner = (Cluster*) DEL_MARK;
                return std::make_pair(parent->parent, parent);
            }
        }

        int fanout = center->get_degree() - parent->get_degree() - children.size();
        if (center == max) fanout++;
        int degree = parent->degree;

        if (fanout < 4 && degree < 4) {
            center->for_all_neighbors([&] (auto neighbor) {
                if(!neighbor->partner && AtomicLoad(&neighbor->parent) == parent) {
                    thread_local_next_root_clusters->push_back(neighbor);
                    AtomicStore(&neighbor->parent, (Cluster*) NULL_PTR);
                }
            });
            if (!center->partner) {
                thread_local_next_root_clusters->push_back(center);
                AtomicStore(&center->parent, (Cluster*) NULL_PTR);
            }
            parent->partner = (Cluster*) DEL_MARK;
        } else {
            if (center == max) {
                if (children.size() > 500) {
                    parlay::parallel_for(0, children.size(), [&] (size_t i) {
                        if (!children[i].second->partner && children[i].second != max) {
                            thread_local_next_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    });
                } else {
                    for (size_t i = 0; i < children.size(); ++i) {
                        if (!children[i].second->partner && children[i].second != max) {
                            thread_local_next_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    }
                }
            } else {
                if (children.size() > 500) {
                    parlay::parallel_for(0, children.size(), [&] (size_t i) {
                        if (!children[i].second->partner) {
                            thread_local_next_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    });
                } else {
                    for (size_t i = 0; i < children.size(); ++i) {
                        if (!children[i].second->partner) {
                            thread_local_next_root_clusters->push_back(children[i].second);
                            AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                        }
                    }
                }
            }
        }

        return std::make_pair(parent->parent, parent);
    });
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_root_clusters(UpdateType update_type) {
    // This function sets the partner fields for all root clusters. For a non-root
    // cluster, we mark its partner field with NON_ROOT_MARK. For root clusters that
    // don't combine with anything, we leave its partner field empty. For high degree
    // root clusters, we assign its partner field as NEW_PAR_MARK, and we add a parent
    // for it in this part. All other root clusters receive no parent at this point.
    // This returns the parent of any non-root clusters that were partnered with.
    thread_local_root_clusters->for_all([&] (Cluster* cluster) {
        if (cluster->get_degree() == 1) {
            recluster_degree_one_root(cluster, update_type);
        }
        else if (cluster->get_degree() == 2) {
            recluster_degree_two_root(cluster);
        }
        else if (cluster->get_degree() >= 3) {
            recluster_high_degree_root(cluster);
        }
    });
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_degree_one_root(Cluster* cluster, UpdateType update_type) {
    // Always partner with a degree 1 or 3+ neighbor. For degree 2
    // neighbors, only attempt to partner with it if it is a non-root
    // cluster that does not already contract. Partnering with degree
    // 2 root clusters will be handled from the deg 2 cluster's side.
    // Return the parent of a non-root combination, if any.
    Cluster* neighbor = cluster->get_neighbor();
    if (neighbor->get_degree() == 1) { // Combine deg 1 root clusters with deg 1 root or non-root clusters
        cluster->partner = neighbor;
        if (neighbor->parent) {
            neighbor->partner = (Cluster*) NON_ROOT_MARK;
            thread_local_del_clusters.push_back(std::make_pair(neighbor->parent->parent, neighbor->parent));
        }
    }
    else if (neighbor->get_degree() == 2) { // Combine deg 1 root cluster with deg 2 non-root clusters that don't contract
        if (neighbor->parent && !neighbor->contracts()) {
            if (CAS(&neighbor->partner, (Cluster*) NULL_PTR, (Cluster*) NON_ROOT_MARK)) {
                cluster->partner = neighbor;
                thread_local_del_clusters.push_back(std::make_pair(neighbor->parent->parent, neighbor->parent));
            }
        }
    }
    else { // Combine deg 1 root cluster with possible deg 3+ non-root clusters
        cluster->partner = neighbor;
        if (update_type == DELETE) {
            if (AtomicLoad(&neighbor->parent))
                if (CAS(&neighbor->partner, (Cluster*) NULL_PTR, (Cluster*) NON_ROOT_MARK))
                    thread_local_del_clusters.push_back(std::make_pair(neighbor->parent->parent, neighbor->parent));
        }
    }
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_degree_two_root(Cluster* cluster) {
    // Only local maxima in priority with respect to deg 2 clusters will act
    if (!is_local_max(cluster)) return;
    // Travel left/right and pair clusters until a deg 3+, deg 1, non-root, or partnered cluster is found
    Cluster* neighbor1 = cluster->get_neighbor();
    Cluster* neighbor2 = cluster->get_other_neighbor(neighbor1);
    for (bool direction : {0, 1}) {
        Cluster* curr = cluster;
        Cluster* next = direction ? neighbor1 : neighbor2;
        if (AtomicLoad(&curr->partner)) {
            curr = next;
            next = curr->get_other_neighbor(cluster);
        }
        while (curr && curr->get_degree() == 2 && !curr->parent && next && next->get_degree() < 3 && !next->contracts()) {
            if (curr != cluster && is_local_max(curr)) break;
            if (next->get_degree() == 2 && !next->parent && is_local_max(next)) break;
            if (!CAS(&curr->partner, (Cluster*) NULL_PTR, next)) break;
            if (next->get_degree() == 1) { // If next deg 1 they can definitely combine
                if (!next->parent) next->partner = curr;
                else next->partner = (Cluster*) NON_ROOT_MARK;
            } else { // deg 2
                Cluster* new_partner = next->parent ? (Cluster*) NON_ROOT_MARK : curr;
                if (!CAS(&next->partner, (Cluster*) NULL_PTR, new_partner)) { // If the CAS fails next was combined from the other side
                    if (AtomicLoad(&next->partner) != curr) // Other side combined with the opposite cluster (you got left hanging)
                        AtomicStore(&curr->partner, (Cluster*) NULL_PTR);
                    break;
                }
            }
            // Both CAS's succeeded or next was degree 1
            if (next->parent) { // Stop traversing at a non-root cluster
                thread_local_del_clusters.push_back(std::make_pair(next->parent->parent, next->parent));
                break;
            }
            if (next->get_degree() == 1) // Stop traversing at deg 1 cluster
                break;
            // Get the next two clusters in the chain
            curr = next->get_other_neighbor(curr);
            if (curr) next = curr->get_other_neighbor(next);
            else break;
        }
    }
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_high_degree_root(Cluster* cluster) {
    // Create the new parent for a high degree root cluster.
    // Find at most one possible non-root degree 1 neighbor,
    // combine with it, and return its parent as del cluster.
    Cluster* parent = allocate_cluster();
    parent->aggregate_weight = cluster->aggregate_weight;
    AtomicStore(&cluster->parent, parent);
    AtomicStore(&cluster->partner, (Cluster*) NEW_PAR_MARK);
    cluster->for_all_neighbors([&] (auto neighbor) {
        if (neighbor->get_degree() == 1 && neighbor->parent) {
            neighbor->parent->partner = (Cluster*) DEL_MARK;
            thread_local_del_clusters.push_back(std::make_pair(neighbor->parent->parent, neighbor->parent));
            neighbor->parent = parent;
        }
    });
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::is_local_max(Cluster* c) {
    // Assumes the input is a degree 2 cluster
    Cluster* neighbor1 = c->get_neighbor();
    Cluster* neighbor2 = c->get_other_neighbor(neighbor1);
    uint64_t hash = hash64((uintptr_t) c);
    uint64_t hash1 = hash64((uintptr_t) neighbor1);
    uint64_t hash2 = hash64((uintptr_t) neighbor2);
    if (neighbor1->get_degree() == 2 && !neighbor1->parent)
        if (hash1 > hash || (hash1 == hash && neighbor1 > c))
            return false;
    if (neighbor2->get_degree() == 2 && !neighbor2->parent)
        if (hash2 > hash || (hash2 == hash && neighbor2 > c))
            return false;
    return true;
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::create_new_parents() {
    // Only returns the brand new clusters to be root clusters at the next level,
    // not the parents of non-root clusters. Those may become root clusters also,
    // but it will be determined by the later step which checks if the grandparent
    // should be deleted.
    if (thread_local_root_clusters->size() > 500) {
        thread_local_root_clusters->for_all([&] (Cluster* cluster) {
            if (cluster->get_degree() == 0) return;
            if (cluster->get_degree() >= 3 && cluster->partner == (Cluster*) NEW_PAR_MARK) {
                thread_local_next_root_clusters->push_back(cluster->parent);
                return;
            }
            Cluster* partner = cluster->partner;
            if (partner) {
                if (partner->partner != cluster) { // Non-root partner or high-degree partner with no partner field set
                    cluster->parent = partner->parent;
                    cluster->parent->aggregate_weight = recompute_parent_aggregate_from_child(cluster);
                } else if (cluster < partner) { // Tie-break for two partnered root clusters
                    Cluster* parent = allocate_cluster();
                    cluster->parent = parent;
                    partner->parent = parent;
                    parent->aggregate_weight = recompute_parent_aggregate_from_pair(cluster, partner);
                    thread_local_next_root_clusters->push_back(parent);
                }
            }
            else { // Non-combining root cluster gets its own parent
                Cluster* parent = allocate_cluster();
                cluster->parent = parent;
                parent->aggregate_weight = cluster->aggregate_weight;
                thread_local_next_root_clusters->push_back(parent);
            }
        });
    } else {
        thread_local_root_clusters->for_all_seq([&] (Cluster* cluster) {
            if (cluster->get_degree() == 0) return;
            if (cluster->get_degree() >= 3 && cluster->partner == (Cluster*) NEW_PAR_MARK) {
                thread_local_next_root_clusters->push_back(cluster->parent);
                return;
            }
            Cluster* partner = cluster->partner;
            if (partner) {
                if (partner->partner != cluster) { // Non-root partner or high-degree partner with no partner field set
                    cluster->parent = partner->parent;
                    cluster->parent->aggregate_weight = recompute_parent_aggregate_from_child(cluster);
                } else if (cluster < partner) {
                    Cluster* parent = allocate_cluster();
                    cluster->parent = parent;
                    partner->parent = parent;
                    parent->aggregate_weight = recompute_parent_aggregate_from_pair(cluster, partner);
                    thread_local_next_root_clusters->push_back(parent);
                }
            }
            else { // Non-combining root cluster gets its own parent
                Cluster* parent = allocate_cluster();
                cluster->parent = parent;
                parent->aggregate_weight = cluster->aggregate_weight;
                thread_local_next_root_clusters->push_back(parent);
            }
        });
    }
}

template<typename aug_t>
void ParallelUFOTree<aug_t>::finish_reclustering() {
    // Clear the partner fields of level i root cluster's and any partnered non-root clusters.
    // For all level i root clusters, insert each incident edge into level i+1 if possible.
    // This code uses the fine-grained locking insert and should work well for only low degree cases.
    if (thread_local_root_clusters->size() > 500) {
        thread_local_root_clusters->for_all([&] (auto cluster) {
            if (!cluster->parent) return; // Only deg 0

            // Clear partner pointers
            Cluster* partner = AtomicLoad(&cluster->partner);
            if (partner == (Cluster*) NEW_PAR_MARK) {
                AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
            }
            else if (partner) {
                Cluster* partner_partner = AtomicLoad(&partner->partner);
                if (partner_partner != cluster) { // Non-root partner
                    if (partner_partner != (Cluster*) NULL_PTR) AtomicStore(&partner->partner, (Cluster*) NULL_PTR);
                    AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
                } else if (cluster < partner) { // Tie-break
                    AtomicStore(&partner->partner, (Cluster*) NULL_PTR);
                    AtomicStore(&cluster->partner, (Cluster*) NULL_PTR);
                }
            }

            // Fill adjacency lists at the next level up
            cluster->for_all_neighbors([&] (auto neighbor) {
                if (neighbor->parent != cluster->parent) {
                    cluster->parent->insert_neighbor(neighbor->parent);
                    neighbor->parent->insert_neighbor(cluster->parent);
                }
            });
        });
    } else {
        thread_local_root_clusters->for_all_seq([&] (auto cluster) {
            if (!cluster->parent) return; // Only deg 0
            // Clear partner pointers
            if (cluster->partner != (Cluster*) NULL_PTR && cluster->partner != (Cluster*) NEW_PAR_MARK)
                cluster->partner->partner = (Cluster*) NULL_PTR;
            cluster->partner = (Cluster*) NULL_PTR;
            // Fill adjacency lists at the next level up
            cluster->for_all_neighbors_seq([&] (auto neighbor) {
                if (neighbor->parent != cluster->parent) {
                    cluster->parent->insert_neighbor(neighbor->parent);
                    neighbor->parent->insert_neighbor(cluster->parent);
                }
            });
        });
    }
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOTree<aug_t>::allocate_cluster() {
    return allocator::create();
    // return new Cluster();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::free_cluster(Cluster* c) {
    c->aggregate_weight = 0;
    allocator::free(c);
    // delete c;
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::UpdateWeight(vertex_t v, uint32_t weight) {
    assert(v >= 0 && v < (vertex_t) leaves.size());
    Cluster* curr = &leaves[v];
    curr->aggregate_weight = weight;
    while (curr->parent) {
        Cluster* parent = curr->parent;
        parent->aggregate_weight = recompute_parent_aggregate_from_child(curr);
        curr = parent;
    }
}

template <typename aug_t>
uint32_t ParallelUFOTree<aug_t>::PathQuery(vertex_t u, vertex_t v) {
    assert(u < leaves.size() && u >= 0 && v < leaves.size() && v >= 0 && connected(u, v));
    if (u == v) return leaves[u].aggregate_weight;

    auto neighbor_at = [] (Cluster* c, int i) -> Cluster* {
        if (!c) return nullptr;
        Cluster* first = c->get_neighbor();
        if (i == 0) return first;
        return c->get_other_neighbor(first);
    };

    uint32_t path_u1, path_u2, path_v1, path_v2;
    path_u1 = path_u2 = leaves[u].aggregate_weight;
    path_v1 = path_v2 = leaves[v].aggregate_weight;
    Cluster *bdry_u1, *bdry_u2, *bdry_v1, *bdry_v2;
    bdry_u1 = bdry_u2 = bdry_v1 = bdry_v2 = nullptr;
    if (leaves[u].get_degree() == 2) {
        bdry_u1 = neighbor_at(&leaves[u], 0);
        bdry_u2 = neighbor_at(&leaves[u], 1);
    }
    if (leaves[v].get_degree() == 2) {
        bdry_v1 = neighbor_at(&leaves[v], 0);
        bdry_v2 = neighbor_at(&leaves[v], 1);
    }
    auto curr_u = &leaves[u];
    auto curr_v = &leaves[v];
    while (curr_u->parent != curr_v->parent) {
        if (curr_u->get_degree() > 2) {
            if (curr_u->parent->get_degree() == 2) {
                // Superunary to Binary
                bdry_u1 = neighbor_at(curr_u->parent, 0);
                bdry_u2 = neighbor_at(curr_u->parent, 1);
                path_u2 = path_u1;
            }
            // Else no changes to representative paths.
        } else {
            for (int i = 0; i < 2; i++) {
                auto neighbor = neighbor_at(curr_u, i);
                if (neighbor && neighbor->parent == curr_u->parent) { // Find the contracting neighbor
                    if (curr_u->get_degree() == 2) {
                        if (curr_u->parent->get_degree() == 2) {
                            // Binary to Binary
                            if (neighbor == bdry_u1) {
                                path_u1 += neighbor->aggregate_weight;
                                bdry_u2 = bdry_u2->parent;
                                for (int i = 0; i < 2; i++)
                                    if (neighbor_at(curr_u->parent, i) && neighbor_at(curr_u->parent, i) != bdry_u2)
                                        bdry_u1 = neighbor_at(curr_u->parent, i);
                            } else {
                                path_u2 += neighbor->aggregate_weight;
                                bdry_u1 = bdry_u1->parent;
                                for (int i = 0; i < 2; i++)
                                    if (neighbor_at(curr_u->parent, i) && neighbor_at(curr_u->parent, i) != bdry_u1)
                                        bdry_u2 = neighbor_at(curr_u->parent, i);
                            }
                        } else {
                            // Binary to Unary
                            path_u1 = (neighbor == bdry_u1) ? path_u2 : path_u1;
                        }
                    } else {
                        if (curr_u->parent->get_degree() == 2) {
                            // Unary to Binary
                            path_u1 = path_u2 = path_u1 + neighbor->aggregate_weight;
                            bdry_u1 = neighbor_at(curr_u->parent, 0);
                            bdry_u2 = neighbor_at(curr_u->parent, 1);
                        } else {
                            // Unary to Unary and Unary to Superunary
                            path_u1 += neighbor->aggregate_weight;
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
        if (curr_v->get_degree() > 2) {
            if (curr_v->parent->get_degree() == 2) {
                // Superunary to Superunary/Binary
                bdry_v1 = neighbor_at(curr_v->parent, 0);
                bdry_v2 = neighbor_at(curr_v->parent, 1);
                path_v2 = path_v1;
            }
            // Else no changes to representative paths.
        } else {
            for (int i = 0; i < 2; i++) {
                auto neighbor = neighbor_at(curr_v, i);
                if (neighbor && neighbor->parent == curr_v->parent) { // Find the contracting neighbor
                    if (curr_v->get_degree() == 2) {
                        if (curr_v->parent->get_degree() == 2) {
                            // Binary to Binary
                            if (neighbor == bdry_v1) {
                                path_v1 += neighbor->aggregate_weight;
                                bdry_v2 = bdry_v2->parent;
                                for (int i = 0; i < 2; i++)
                                    if (neighbor_at(curr_v->parent, i) && neighbor_at(curr_v->parent, i) != bdry_v2)
                                        bdry_v1 = neighbor_at(curr_v->parent, i);
                            } else {
                                path_v2 += neighbor->aggregate_weight;
                                bdry_v1 = bdry_v1->parent;
                                for (int i = 0; i < 2; i++)
                                    if (neighbor_at(curr_v->parent, i) && neighbor_at(curr_v->parent, i) != bdry_v1)
                                        bdry_v2 = neighbor_at(curr_v->parent, i);
                            }
                        } else {
                            // Binary to Unary
                            path_v1 = (neighbor == bdry_v1) ? path_v2 : path_v1;
                        }
                    } else {
                        if (curr_v->parent->get_degree() == 2) {
                            // Unary to Binary
                            path_v1 = path_v2 = path_v1 + neighbor->aggregate_weight;
                            bdry_v1 = neighbor_at(curr_v->parent, 0);
                            bdry_v2 = neighbor_at(curr_v->parent, 1);
                        } else {
                            // Unary to Unary and Unary to Superunary
                            path_v1 += neighbor->aggregate_weight;
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
    uint32_t total = 0;
    if (curr_u->get_degree() == 2)
        total += curr_v == bdry_u1 ? path_u1 : path_u2;
    else
        total += path_u1;
    if (curr_v->get_degree() == 2)
        total += curr_u == bdry_v1 ? path_v1 : path_v2;
    else
        total += path_v1;
    // If the LCA contracts them in a star mergeadd the weight of the center vertex
    if (curr_u->get_degree() == 1 && curr_v->get_degree() == 1
    && neighbor_at(curr_u, 0) != curr_v) [[unlikely]] {
        total += neighbor_at(curr_u, 0)->aggregate_weight;
    }
    return total;
}

// template <typename aug_t>
// uint32_t ParallelUFOTree<aug_t>::PathQuery(vertex_t u, vertex_t v) {
//     assert(u < leaves.size() && u >= 0 && v < leaves.size() && v >= 0 && connected(u, v));
//     if (u == v) return leaves[u].aggregate_weight;

//     uint32_t path_u1 = 0, path_u2 = 0, path_v1 = 0, path_v2 = 0;
//     Cluster *bdry_u1 = nullptr, *bdry_u2 = nullptr, *bdry_v1 = nullptr, *bdry_v2 = nullptr;

//     if (leaves[u].get_degree() == 2) {
//         bdry_u1 = leaves[u].get_neighbor();
//         bdry_u2 = leaves[u].get_other_neighbor(bdry_u1);
//     }
//     if (leaves[v].get_degree() == 2) {
//         bdry_v1 = leaves[v].get_neighbor();
//         bdry_v2 = leaves[v].get_other_neighbor(bdry_v1);
//     }

//     Cluster* curr_u = &leaves[u];
//     Cluster* curr_v = &leaves[v];
//     while (curr_u->parent != curr_v->parent) {
//         if (curr_u->get_degree() > 2) {
//             if (curr_u != &leaves[u] && curr_u != &leaves[v]
//                 && !is_endpoint_subtree_aggregate(leaves, curr_u, u, v))
//                 path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//             if (curr_u->parent->get_degree() == 2) {
//                 bdry_u1 = curr_u->parent->get_neighbor();
//                 bdry_u2 = curr_u->parent->get_other_neighbor(bdry_u1);
//                 path_u2 = path_u1;
//             }
//         } else {
//             Cluster* n1 = curr_u->get_neighbor();
//             Cluster* n2 = curr_u->get_other_neighbor(n1);
//             bool contracted_u = false;
//             for (Cluster* neighbor : {n1, n2}) {
//                 if (neighbor && neighbor->parent == curr_u->parent) {
//                     contracted_u = true;
//                     if (curr_u->get_degree() == 2) {
//                         if (curr_u->parent->get_degree() == 2) {
//                             if (neighbor == bdry_u1) {
//                                 path_u1 = combine_weights(path_u1, neighbor->aggregate_weight);
//                                 if (curr_u != &leaves[u] && curr_u != &leaves[v] && curr_u->aggregate_weight > 0)
//                                     path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//                                 bdry_u2 = bdry_u2 ? bdry_u2->parent : nullptr;
//                                 Cluster* p1 = curr_u->parent->get_neighbor();
//                                 Cluster* p2 = curr_u->parent->get_other_neighbor(p1);
//                                 bdry_u1 = (p1 && p1 != bdry_u2) ? p1 : p2;
//                             } else {
//                                 path_u2 = combine_weights(path_u2, neighbor->aggregate_weight);
//                                 if (curr_u != &leaves[u] && curr_u != &leaves[v] && curr_u->aggregate_weight > 0)
//                                     path_u2 = combine_weights(path_u2, curr_u->aggregate_weight);
//                                 bdry_u1 = bdry_u1 ? bdry_u1->parent : nullptr;
//                                 Cluster* p1 = curr_u->parent->get_neighbor();
//                                 Cluster* p2 = curr_u->parent->get_other_neighbor(p1);
//                                 bdry_u2 = (p1 && p1 != bdry_u1) ? p1 : p2;
//                             }
//                         } else {
//                             path_u1 = (neighbor == bdry_u1) ? path_u1 : path_u2;
//                             if (curr_u->get_degree() == 2 && curr_u->parent && curr_u->parent->get_degree() == 1
//                                 && curr_u != &leaves[u] && curr_u != &leaves[v] && curr_u->aggregate_weight > 0
//                                 && curr_u->aggregate_weight != leaves[u].aggregate_weight
//                                 && curr_u->aggregate_weight != leaves[v].aggregate_weight
//                                 && !is_endpoint_leaf_parent_duplicate(leaves, curr_u, u, v)
//                                 && !is_endpoint_subtree_aggregate(leaves, curr_u, u, v))
//                                 path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//                         }
//                     } else {
//                         if (curr_u->parent->get_degree() == 2) {
//                             if (curr_u != &leaves[u])
//                                 path_u1 = path_u2 = combine_weights(path_u1, curr_u->aggregate_weight);
//                             bdry_u1 = curr_u->parent->get_neighbor();
//                             bdry_u2 = curr_u->parent->get_other_neighbor(bdry_u1);
//                         } else if (curr_u->parent->get_degree() > 2) {
//                             if (curr_u != &leaves[u])
//                                 path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//                         } else {
//                             path_u1 = combine_weights(path_u1, neighbor->aggregate_weight);
//                         }
//                     }
//                     break;
//                 }
//             }
//             if (!contracted_u && curr_u != &leaves[u] && curr_u != &leaves[v] && curr_u->aggregate_weight > 0
//                 && curr_u->aggregate_weight != leaves[u].aggregate_weight
//                 && curr_u->aggregate_weight != leaves[v].aggregate_weight
//                 && !is_endpoint_leaf_parent_duplicate(leaves, curr_u, u, v)
//                 && !is_endpoint_subtree_aggregate(leaves, curr_u, u, v)) {
//                 if (curr_u->get_degree() == 2 && curr_u->parent && curr_u->parent->get_degree() == 1)
//                     path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//                 else if (curr_u->get_degree() == 1 && curr_u->parent && curr_u->parent->get_degree() == 2)
//                     path_u1 = combine_weights(path_u1, curr_u->aggregate_weight);
//             }
//             if (!curr_u->contracts()) {
//                 if (bdry_u1) bdry_u1 = bdry_u1->parent;
//                 if (bdry_u2) bdry_u2 = bdry_u2->parent;
//             }
//         }
//         curr_u = curr_u->parent;

//         if (curr_v->get_degree() > 2) {
//             if (curr_v != &leaves[v] && curr_v != &leaves[u]
//                 && !is_endpoint_subtree_aggregate(leaves, curr_v, u, v))
//                 path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//             if (curr_v->parent->get_degree() == 2) {
//                 bdry_v1 = curr_v->parent->get_neighbor();
//                 bdry_v2 = curr_v->parent->get_other_neighbor(bdry_v1);
//                 path_v2 = path_v1;
//             }
//         } else {
//             Cluster* n1 = curr_v->get_neighbor();
//             Cluster* n2 = curr_v->get_other_neighbor(n1);
//             bool contracted_v = false;
//             for (Cluster* neighbor : {n1, n2}) {
//                 if (neighbor && neighbor->parent == curr_v->parent) {
//                     contracted_v = true;
//                     if (curr_v->get_degree() == 2) {
//                         if (curr_v->parent->get_degree() == 2) {
//                             if (neighbor == bdry_v1) {
//                                 path_v1 = combine_weights(path_v1, neighbor->aggregate_weight);
//                                 if (curr_v != &leaves[u] && curr_v != &leaves[v] && curr_v->aggregate_weight > 0)
//                                     path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//                                 bdry_v2 = bdry_v2 ? bdry_v2->parent : nullptr;
//                                 Cluster* p1 = curr_v->parent->get_neighbor();
//                                 Cluster* p2 = curr_v->parent->get_other_neighbor(p1);
//                                 bdry_v1 = (p1 && p1 != bdry_v2) ? p1 : p2;
//                             } else {
//                                 path_v2 = combine_weights(path_v2, neighbor->aggregate_weight);
//                                 if (curr_v != &leaves[u] && curr_v != &leaves[v] && curr_v->aggregate_weight > 0)
//                                     path_v2 = combine_weights(path_v2, curr_v->aggregate_weight);
//                                 bdry_v1 = bdry_v1 ? bdry_v1->parent : nullptr;
//                                 Cluster* p1 = curr_v->parent->get_neighbor();
//                                 Cluster* p2 = curr_v->parent->get_other_neighbor(p1);
//                                 bdry_v2 = (p1 && p1 != bdry_v1) ? p1 : p2;
//                             }
//                         } else {
//                             path_v1 = (neighbor == bdry_v1) ? path_v1 : path_v2;
//                             if (curr_v->get_degree() == 2 && curr_v->parent && curr_v->parent->get_degree() == 1
//                                 && curr_v != &leaves[u] && curr_v != &leaves[v] && curr_v->aggregate_weight > 0
//                                 && curr_v->aggregate_weight != leaves[u].aggregate_weight
//                                 && curr_v->aggregate_weight != leaves[v].aggregate_weight
//                                 && !is_endpoint_leaf_parent_duplicate(leaves, curr_v, u, v)
//                                 && !is_endpoint_subtree_aggregate(leaves, curr_v, u, v))
//                                 path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//                         }
//                     } else {
//                         if (curr_v->parent->get_degree() == 2) {
//                             if (curr_v != &leaves[v])
//                                 path_v1 = path_v2 = combine_weights(path_v1, curr_v->aggregate_weight);
//                             bdry_v1 = curr_v->parent->get_neighbor();
//                             bdry_v2 = curr_v->parent->get_other_neighbor(bdry_v1);
//                         } else if (curr_v->parent->get_degree() > 2) {
//                             if (curr_v != &leaves[v])
//                                 path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//                         } else {
//                             path_v1 = combine_weights(path_v1, neighbor->aggregate_weight);
//                         }
//                     }
//                     break;
//                 }
//             }
//             if (!contracted_v && curr_v != &leaves[u] && curr_v != &leaves[v] && curr_v->aggregate_weight > 0
//                 && curr_v->aggregate_weight != leaves[u].aggregate_weight
//                 && curr_v->aggregate_weight != leaves[v].aggregate_weight
//                 && !is_endpoint_leaf_parent_duplicate(leaves, curr_v, u, v)
//                 && !is_endpoint_subtree_aggregate(leaves, curr_v, u, v)) {
//                 if (curr_v->get_degree() == 2 && curr_v->parent && curr_v->parent->get_degree() == 1)
//                     path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//                 else if (curr_v->get_degree() == 1 && curr_v->parent && curr_v->parent->get_degree() == 2)
//                     path_v1 = combine_weights(path_v1, curr_v->aggregate_weight);
//             }
//             if (!curr_v->contracts()) {
//                 if (bdry_v1) bdry_v1 = bdry_v1->parent;
//                 if (bdry_v2) bdry_v2 = bdry_v2->parent;
//             }
//         }
//         curr_v = curr_v->parent;
//     }

//     uint32_t total = 0;
//     total = combine_weights(total, (curr_u->get_degree() == 2) ? ((curr_v == bdry_u1) ? path_u1 : path_u2) : path_u1);
//     total = combine_weights(total, (curr_v->get_degree() == 2) ? ((curr_u == bdry_v1) ? path_v1 : path_v2) : path_v1);

//     if (curr_u == curr_v) {
//         total = combine_weights(total, curr_u->aggregate_weight);
//     } else if (curr_u->get_degree() == 1 && curr_v->get_degree() == 1) {
//         Cluster* nu = curr_u->get_neighbor();
//         Cluster* nv = curr_v->get_neighbor();
//         if (nu != curr_v) {
//             Cluster* center = advance_to_nonzero_cluster(nu);
//             Cluster* other_center = advance_to_nonzero_cluster(nv);
//             if (center && center == other_center) {
//                 total = combine_weights(total, center->aggregate_weight);
//             } else if (curr_u->parent && curr_u->parent == curr_v->parent) {
//                 total = combine_weights(total, curr_u->parent->aggregate_weight);
//             }
//         } else if (curr_u->parent && curr_u->parent == curr_v->parent) {
//             total = combine_weights(total, curr_u->parent->aggregate_weight);
//         }
//     } else {
//         bool final_is_leaf_pair =
//             ((curr_u == &leaves[u] && curr_v == &leaves[v]) ||
//              (curr_u == &leaves[v] && curr_v == &leaves[u]));
//         uint32_t pair_value = final_is_leaf_pair ? 0 : recompute_parent_aggregate_from_pair(curr_u, curr_v);
//         bool query_has_high_degree_endpoint = (leaves[u].get_degree() > 2 || leaves[v].get_degree() > 2);
//         if (pair_value == 0 && query_has_high_degree_endpoint) {
//             if (curr_u->get_degree() == 1 && curr_v->get_degree() == 2 && curr_v != &leaves[u] && curr_v != &leaves[v])
//                 pair_value = curr_v->aggregate_weight;
//             else if (curr_v->get_degree() == 1 && curr_u->get_degree() == 2 && curr_u != &leaves[u] && curr_u != &leaves[v])
//                 pair_value = curr_u->aggregate_weight;
//         }
//         total = combine_weights(total, pair_value);
//     }

//     bool adjacent_deg2_leaf_pair = false;
//     if (curr_u == &leaves[u] && curr_v == &leaves[v] && curr_u->get_degree() == 2 && curr_v->get_degree() == 2) {
//         Cluster* n = curr_u->get_neighbor();
//         Cluster* on = curr_u->get_other_neighbor(n);
//         adjacent_deg2_leaf_pair = (n == curr_v || on == curr_v);
//     }
//     if (!adjacent_deg2_leaf_pair) {
//         total = combine_weights(total, leaves[u].aggregate_weight);
//         total = combine_weights(total, leaves[v].aggregate_weight);
//     }
//     return total;
// }

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}

