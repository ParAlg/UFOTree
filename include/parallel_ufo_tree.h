#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include "aggregator.h"
#include "parallel_ufo_cluster.h"


extern parlay::internal::timer timer0;
extern parlay::internal::timer timer1;
extern parlay::internal::timer timer2;
extern parlay::internal::timer timer3;
extern parlay::internal::timer timer4;
extern parlay::internal::timer timer5;

extern parlay::internal::timer subtimer1;
extern parlay::internal::timer subtimer2;
extern parlay::internal::timer subtimer3;
extern parlay::internal::timer subtimer4;

namespace dgbs {

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
public:
    // UFO tree interface
    ParallelUFOTree(vertex_t n, vertex_t k);
    ~ParallelUFOTree();
    void batch_link(parlay::sequence<std::pair<int, int>>& links);
    void batch_cut(parlay::sequence<std::pair<int, int>>& cuts);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
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
    thread_local_del_edges(), thread_local_new_del_edges() {}

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
    timer0.start();

    parlay::sequence<std::pair<Cluster*, Cluster*>> del_clusters;
    parlay::sequence<std::pair<Cluster*, Cluster*>> next_del_clusters;
    parlay::sequence<parlay::sequence<std::pair<Cluster*, Cluster*>>> all_del_clusters;

    timer1.start();
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
    auto dir_update_groups = integer_group_by_key_inplace(dir_updates);

    parlay::sequence<std::pair<Cluster*, Cluster*>> parents;
    parlay::sequence<std::pair<Cluster*, EdgeSlice>> parent_groups;
    parlay::parallel_do(
    [&] () {
        // Group the affected vertices by parent.
        parents = parlay::map(dir_update_groups, [&] (auto group) {
            return std::make_pair(group.first->parent, group.first);
        });
        parent_groups = integer_group_by_key_inplace(parents);
    },
    [&] () {
        // For deletion batches, keep a trace of the level i+2 del edges to delete and decrement `degree` at each level.
        if (update_type == DELETE) {
            auto del_edges = thread_local_del_edges.to_sequence();
            thread_local_del_edges.clear();
            auto del_edge_groups = integer_group_by_key_inplace(del_edges);
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
    timer1.stop();


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

        timer2.start();
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
        dir_update_groups = integer_group_by_key_inplace(dir_updates);
        timer2.stop();

        // =======
        // PHASE 2
        // =======

        timer3.start();
        // Recluster the root clusters.
        recluster_root_clusters(update_type);
        del_clusters = parlay::append(del_clusters, thread_local_del_clusters.to_sequence());
        thread_local_del_clusters.clear();

        // This returns only the new clusters that were created during the reclustering at this level.
        create_new_parents();
        timer3.stop();

        // =======
        // PHASE 3
        // =======

        timer4.start();
        subtimer1.start();
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
        subtimer1.stop();

        subtimer2.start();
        parlay::parallel_do(
        [&] () {
            // Group edges to delete by endpoint.
            new_del_edge_groups = integer_group_by_key_inplace(new_del_edges);
        },
        [&] () {
            parlay::parallel_do(
            [&] () {
                // Group del clusters by parent.
                parent_groups = integer_group_by_key_inplace(del_clusters);
            },
            [&] () {
                // For deletion batches, keep a trace of the level i+2 del edges to delete and decrement `degree` at each level.
                if (update_type == DELETE) {
                    // Group the level i+2 updates by endpoint.
                    auto del_edges = thread_local_del_edges.to_sequence();
                    thread_local_del_edges.clear();
                    auto del_edge_groups = integer_group_by_key_inplace(del_edges);
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
        subtimer2.stop();

        // The next level i+1 root clusters, which are children of a level i+2 del cluster that will be deleted.
        // This function populates `thread_local_next_root_clusters` with the next level i+1 root clusters.
        // This function returns the next level i+2 del clusters and we store it in `next_del_clusters`.
        subtimer3.start();
        next_del_clusters = process_del_clusters(parent_groups);
        subtimer3.stop();

        // Delete pointers to level i+1 del clusters that will get deleted.
        subtimer4.start();
        parlay::parallel_for(0, new_del_edge_groups.size(), [&] (size_t i) {
            auto& [cluster, edges] = new_del_edge_groups[i];
            auto neighbors = parlay::map(edges, [&] (auto x) { return x.second; });
            cluster->delete_neighbors_sorted(neighbors);
        });
        subtimer4.stop();

        // Delete level i+1 del clusters that should be deleted.
        all_del_clusters.push_back(std::move(del_clusters));
        timer4.stop();

        // =======
        // PHASE 4
        // =======

        // Clear the partner fields and populate neighbors of level i+1 clusters.
        timer5.start();
        finish_reclustering();
        timer5.stop();

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
    timer0.stop();
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

        Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
        int max_degree = max->get_degree();
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
                parlay::parallel_for(0, children.size(), [&] (size_t i) {
                    if (children[i].second != max) {
                        thread_local_root_clusters->push_back(children[i].second);
                        AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                    }
                });
            } else {
                parlay::parallel_for(0, children.size(), [&] (size_t i) {
                    thread_local_root_clusters->push_back(children[i].second);
                    AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                });
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

        Cluster* max = (*parlay::max_element(children, [&] (auto x, auto y) { return x.second->get_degree() < y.second->get_degree(); })).second;
        int max_degree = max->get_degree();
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
                parlay::parallel_for(0, children.size(), [&] (size_t i) {
                    if (!children[i].second->partner && children[i].second != max) {
                        thread_local_next_root_clusters->push_back(children[i].second);
                        AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                    }
                });
            } else {
                parlay::parallel_for(0, children.size(), [&] (size_t i) {
                    if (!children[i].second->partner) {
                        thread_local_next_root_clusters->push_back(children[i].second);
                        AtomicStore(&children[i].second->parent, (Cluster*) NULL_PTR);
                    }
                });
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
            } else if (cluster < partner) { // Tie-break for two partnered root clusters
                Cluster* parent = allocate_cluster();
                cluster->parent = parent;
                partner->parent = parent;
                thread_local_next_root_clusters->push_back(parent);
            }
        }
        else { // Non-combining root cluster gets its own parent
            Cluster* parent = allocate_cluster();
            cluster->parent = parent;
            thread_local_next_root_clusters->push_back(parent);
        }
    });
}

template<typename aug_t>
void ParallelUFOTree<aug_t>::finish_reclustering() {
    // Clear the partner fields of level i root cluster's and any partnered non-root clusters.
    // For all level i root clusters, insert each incident edge into level i+1 if possible.
    // This code uses the fine-grained locking insert and should work well for only low degree cases.
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
}

template <typename aug_t>
ParallelUFOCluster<aug_t>* ParallelUFOTree<aug_t>::allocate_cluster() {
    return allocator::create();
    // return new Cluster();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::free_cluster(Cluster* c) {
    allocator::free(c);
    // delete c;
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}
