#pragma once
#include <parlay/parallel.h>
#include <parlay/sequence.h>
#include "types.h"
#include "util.h"
#include "bridge.h"
#include "hash_bag.h"
#include "parallel_ufo_cluster.h"


namespace dgbs {

template <typename aug_t = empty_t>
class ParallelUFOTree {
    using Cluster = ParallelUFOCluster<aug_t>;
    using allocator = parlay::type_allocator<Cluster>;
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
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    // Helper functions
    void recluster_tree(parlay::sequence<std::pair<int, int>>& updates, UpdateType update_type);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : leaves(n) {}

template <typename aug_t>
ParallelUFOTree<aug_t>::~ParallelUFOTree() {
    allocator::finish();
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
    parlay::sequence<Cluster*> root_clusters;
    parlay::sequence<Cluster*> del_clusters;
    parlay::sequence<std::pair<Cluster*, Cluster*>> dir_updates;

    // The intial dir updates are just the batch of links or cuts in both directions.
    dir_updates = parlay::tabulate(2*updates.size(), [&] (size_t i) {
        if (i % 2 == 0) return std::make_pair(&leaves[updates[i/2].first], &leaves[updates[i/2].second]);
        return std::make_pair(&leaves[updates[i/2].second], &leaves[updates[i/2].first]);
    });
    
    // Group by endpoint. Then group by parent.
    auto dir_update_groups = parlay::group_by_key(dir_updates);
    auto parents = parlay::delayed_tabulate(dir_update_groups.size(), [&] (size_t i) {
        return std::make_pair(dir_update_groups[i].first->parent, dir_update_groups[i].first);
    });
    auto parent_groups = parlay::group_by_key(parents);

    // The level 0 root clusters are children of level 1 parents that will get deleted, or deg 1 endpoints.
    root_clusters = parlay::flatten(parlay::delayed_tabulate(parent_groups.size(), [&] (size_t i) {
        auto& [parent, children] = parent_groups[i];
        parlay::sequence<Cluster*> local_root_clusters;
        if (parent == nullptr) return children;

        Cluster* max = *parlay::max_element(children, [&] (Cluster* x, Cluster* y) { return x->get_degree() < y->get_degree(); });
        size_t max_degree = max->get_degree();

        if (max_degree == 1) {
            Cluster* center = max->get_neighbor();
            if (center->parent != max->parent) {
                local_root_clusters.push_back(max);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            } else {
                size_t adjusted_center_degree = center->get_degree() - children.size();
                local_root_clusters = std::move(children);
                if (adjusted_center_degree < 3) {
                    center->delete_neighbors(local_root_clusters); // Should replace delete / reinsert with parallel iteration
                    if (adjusted_center_degree >= 1) {
                        Cluster* neighbor1 = center->get_neighbor();
                        if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                        if (adjusted_center_degree >= 2) {
                            Cluster* neighbor2 = center->get_other_neighbor(neighbor1);
                            if (neighbor2->parent == parent) local_root_clusters.push_back(center->get_other_neighbor(neighbor2));
                        }
                    }
                    center->insert_neighbors(local_root_clusters);
                    local_root_clusters.push_back(center);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                }
            }
        }

        else if (max_degree == 2) {
            local_root_clusters.push_back(max);
            Cluster* neighbor1 = max->get_neighbor();
            Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
            if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
            else local_root_clusters.push_back(neighbor2);
            parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
        }

        else if (max_degree >= 3) {
            local_root_clusters = std::move(children);
            if (max_degree - (children.size()-1) < 3) {
                max->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                Cluster* neighbor1 = max->get_neighbor();
                if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                if (max_degree - children.size() == 2) {
                    Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                    if (neighbor2->parent == parent) local_root_clusters.push_back(max->get_other_neighbor(neighbor2));
                }
                max->insert_neighbors(children);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            } else {
                local_root_clusters.erase(parlay::find(local_root_clusters, max));
            }
        }

        parlay::parallel_for(0, local_root_clusters.size(), [&] (size_t i) {
            local_root_clusters[i]->parent = nullptr;
        });
        return local_root_clusters;
    }));

    // The intial del clusters are just all parents of vertices that got an update.
    del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
        if (!group.first) return std::nullopt;
        return group.first;
    });

    /* This is the main loop over all levels of the tree. At each level i, we assume that `root_clusters`is
    populated with the level i clusters formed from contractions in the previous level, and that `del_clusters`
    is populated with the level i+1 clusters that were parents of del clusters at the previous level, along
    with the parents of any nodes that started a new remove ancestor path due to a root cluster merging into a
    non-root cluster that didn't already contract.
    (1) First we take the level i+1 `del_clusters`, and populate `next_del_clusters` by mapping them to their
    level i+2 parents. We also determine which ones will actually be deleted. During this step we delete the
    del clusters. We also find any level i+1 non-del children of the level i+2 clusters which will be deleted.
    These will become `next_root_clusters` at the next level i+1.
    (2) For del clusters that don't get deleted, we need to update their adjacency lists based on the batch of
    links or cuts. For cuts it is necessary to do this after determining the `next_root_clusters` in step 1, so
    that the siblings of a cluster can still be found.
    (3) Next we, will determine the new set of contractions over the level i root clusters using an MIS. While
    doing this we also allocate new level i+1 parent nodes and update the parent pointers of the root clusters.
    These new parents become the `next_root_clusters` along with those determined in step 1. Any contraction
    with a non-root cluster adds its level i+2 grandparent to `next_del_clusters`. After, we fill the adjacency
    lists of our new parent clusters and add them to the adjacency list of their neighbors. */
    while (root_clusters.size() > 0 || del_clusters.size() > 0) {
        // =======
        // PHASE 2
        // =======

        // Insert or delete the edges at this level.
        auto dir_update_groups = parlay::group_by_key(dir_updates);
        parlay::parallel_for(0, dir_update_groups.size(), [&] (size_t i) {
            auto& [cluster, neighbors] = dir_update_groups[i];
            if (update_type == INSERT) cluster->insert_neighbors(neighbors);
            else cluster->delete_neighbors(neighbors);
        });

        // Map the dir updates to the next level.
        dir_updates = parlay::map_maybe(dir_updates, [&] (std::pair<Cluster*, Cluster*> x) -> std::optional<std::pair<Cluster*, Cluster*>> {
            if (x.first->parent != x.second->parent)
                return std::make_pair(x.first->parent, x.second->parent);
            return std::nullopt;
        });

        // =======
        // PHASE 3
        // =======

        /* We also need to handle the case where a root cluster became high degree because of many
        insertions, but it needs to find a degree 1 neighbor that is a non-root cluster because it
        previously happened to not contract with this cluster. Fortunately we can detect when this
        will happen by looking at the number of edges inserted, and we can grab that deg 1 neighbor
        easily before inserting all the edges and making this cluster high degree.
        This means populating the adjacency lists should happen during phase 3 of the next level. */


        // We can probably eliminate the extra `partner` field per cluster by using the `parent` field smartly.
        // This returns the parent's of non-root clusters that were combined with, starting new RA paths.
        auto additional_del_clusters = parlay::flatten(parlay::delayed_tabulate(root_clusters.size(), [&] (size_t i) {
            Cluster* cluster = root_clusters[i];
            parlay::sequence<Cluster*> local_del_clusters;
            if (cluster->get_degree() == 1) {
                Cluster* neighbor = cluster->get_neighbor();
                // Combine deg 1 root clusters with deg 1 root or non-root clusters
                if (neighbor->get_degree() == 1) {
                    cluster->partner = neighbor;
                    if (neighbor->parent && !neighbor->partner) {
                        neighbor->partner = neighbor; // Mark non-root contracting clusters with a self-partner
                        cluster->parent = neighbor->parent;
                        local_del_clusters.push_back(neighbor->parent);
                    }
                    else if (!neighbor->parent && cluster < neighbor) {
                        Cluster* parent = allocator::create();
                        cluster->parent = parent;
                        neighbor->parent = parent;
                    }
                }
                // Combine deg 1 root cluster with deg 2 non-root clusters that don't combine
                else if (neighbor->get_degree() == 2) {
                    if (neighbor->parent && !neighbor->partner) {
                        if (!neighbor->contracts()) {
                            // Mark non-root contracting clusters with a self-partner
                            if (!neighbor->partner && gbbs::CAS(&neighbor->partner, (Cluster*) nullptr, neighbor)) {
                                cluster->partner = neighbor;
                                cluster->parent = neighbor->parent;
                                local_del_clusters.push_back(neighbor->parent);
                            }
                        } else { // Deg 1 cluster neighboring a contracting deg 2 non-root cluster gets its own parent
                            cluster->partner = cluster;
                            Cluster* parent = allocator::create();
                            cluster->parent = parent;
                        }
                    }
                }
                // Combine deg 1 root cluster with deg 3+ clusters always
                else if (neighbor->get_degree() >= 3) {
                    cluster->partner = neighbor;
                    neighbor->partner = cluster;
                    if (!neighbor->parent) {
                        Cluster* parent = allocator::create();
                        if (!gbbs::CAS(&neighbor->parent, (Cluster*) nullptr, parent))
                            allocator::destroy(parent);
                    }
                    else if (neighbor->parent && !neighbor->partner) {
                        local_del_clusters.push_back(neighbor->parent);
                    }
                    cluster->parent = neighbor->parent;
                }
            }
            else if (cluster->get_degree() == 2) {
                // Only local maxima in priority with respect to deg 2 clusters will act
                Cluster* neighbor1 = cluster->get_neighbor();
                Cluster* neighbor2 = cluster->get_other_neighbor(neighbor1);
                uint64_t hash = hash64((uintptr_t) cluster);
                uint64_t hash1 = hash64((uintptr_t) neighbor1);
                uint64_t hash2 = hash64((uintptr_t) neighbor2);
                if (neighbor1->get_degree() == 2 && (!neighbor1->parent || neighbor1->partner))
                    if (hash1 > hash || (hash1 == hash && neighbor1 > cluster))
                        return local_del_clusters;
                if (neighbor2->get_degree() == 2 && (!neighbor2->parent || neighbor2->partner))
                    if (hash2 > hash || (hash2 == hash && neighbor2 > cluster))
                        return local_del_clusters;
                // Travel left/right and pair clusters until a deg 3, deg 1, non-root, or partnered cluster is found
                for (auto direction : {0, 1}) {
                    auto curr = cluster;
                    auto next = direction ? neighbor1 : neighbor2;
                    if (curr->partner) {
                        curr = next;
                        next = curr->get_other_neighbor(cluster);
                    }
                    while (curr && !curr->parent && curr->get_degree() == 2 && next && next->get_degree() < 3 && !next->contracts()) {
                        if (curr->partner || !gbbs::CAS(&curr->partner, (Cluster*) nullptr, next)) break;
                        if (next->get_degree() == 1) { // If next deg 1 they can definitely combine
                            if (!next->parent) next->partner = curr;
                            else next->partner = next; // Mark non-root contracting clusters with a self-partner
                        } else {
                            Cluster* new_partner = next->parent ? next : curr; // Mark non-root contracting clusters with a self-partner
                            if (next->partner || !gbbs::CAS(&next->partner, (Cluster*) nullptr, new_partner)) { // If the CAS fails next was combined from the other side
                                if (next->partner != curr) curr->partner = nullptr;
                                break;
                            }
                        }
                        // Both CAS's succeeded or next was degree 1
                        if (next->parent) {
                            curr->parent = next->parent;
                            local_del_clusters.push_back(next->parent);
                            break; // Stop traversing at a non-root cluster
                        }
                        else {
                            Cluster* parent = allocator::create();
                            curr->parent = parent;
                            next->parent = parent;
                        }
                        if (next->get_degree() == 1) break;
                        // Get the next two clusters in the chain
                        curr = next->get_other_neighbor(curr);
                        if (curr) next = curr->get_other_neighbor(next);
                    }
                    if (curr && !curr->parent && curr->get_degree() == 1) { // Deg 1 cluster neighboring a contracting deg 2 root cluster gets its own parent
                        curr->partner = curr;
                        Cluster* parent = allocator::create();
                        curr->parent = parent;
                    }
                }
                // Deg 2 root cluster that doesn't combine gets its own parent
                if (!cluster->partner) {
                    cluster->partner = cluster;
                    Cluster* parent = allocator::create();
                    cluster->parent = parent;
                }
            }
            return local_del_clusters;
        }));

        // Fill the adjacency lists of new clusters
        // This should work well for only linked list cases
        // This returns only the new clusters that were created during the reclustering at this level.
        auto next_root_clusters1 = parlay::flatten(parlay::delayed_tabulate(root_clusters.size(), [&] (size_t i) {
            auto cluster = root_clusters[i];
            parlay::sequence<Cluster*> local_root_clusters;
            if (!cluster->parent) return local_root_clusters; // Only deg 0

            // Clear partner pointers and find the newly created clusters
            Cluster* partner = cluster->partner;
            if (partner && partner != cluster) {
                if (partner->partner != cluster) { // Non-root partner
                    partner->partner = nullptr;
                    cluster->partner = nullptr;
                } else if (cluster < partner) { // Tie-break
                    cluster->partner->partner = nullptr;
                    cluster->partner = nullptr;
                    local_root_clusters.push_back(cluster->parent);
                }
            } else if (partner) { // Non-combining cluster has its own parent
                cluster->partner = nullptr;
                local_root_clusters.push_back(cluster->parent);
            }

            // Fill adjacency lists at the next level up
            Cluster* neighbor1 = cluster->get_neighbor();
            Cluster* neighbor2 = cluster->get_other_neighbor(neighbor1);
            if (neighbor1 && cluster->parent != neighbor1->parent) {
                cluster->parent->insert_neighbor(neighbor1->parent);
                neighbor1->parent->insert_neighbor(cluster->parent);
            }
            if (neighbor2 && cluster->parent != neighbor2->parent) {
                cluster->parent->insert_neighbor(neighbor2->parent);
                neighbor2->parent->insert_neighbor(cluster->parent);
            }

            return local_root_clusters;
        }));

        // =======
        // PHASE 1
        // =======

        del_clusters.append(additional_del_clusters);

        // Group del clusters by parent.
        auto parents = parlay::delayed_tabulate(del_clusters.size(), [&] (size_t i) {
            return std::make_pair(del_clusters[i]->parent, del_clusters[i]);
        });
        auto parent_groups = parlay::group_by_key(parents);
        // Get the next del clusters.
        auto next_del_clusters = parlay::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
            if (!group.first) return std::nullopt;
            return group.first;
        });
        // Get next root clusters, clear their parent pointers, delete some current del clusters.
        // This should return any del cluster that doesn't get deleted, but whose parent get's deleted.
        // This should also return any cluster that is not a del cluster, but is a child of a deleted parent.
        auto next_root_clusters2 = parlay::flatten(parlay::delayed_tabulate(parent_groups.size(), [&] (size_t i) {
            auto& [parent, children] = parent_groups[i];
            parlay::sequence<Cluster*> local_root_clusters;
            if (parent == nullptr) return local_root_clusters;

            Cluster* max = *parlay::max_element(children, [&] (Cluster* x, Cluster* y) { return x->get_degree() < y->get_degree(); });
            size_t max_degree = max->get_degree();

            if (max_degree == 1) {
                Cluster* center = max->get_neighbor();
                if (center->parent != max->parent) {
                    local_root_clusters.push_back(max);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                } else {
                    size_t adjusted_center_degree = center->get_degree() - children.size();
                    local_root_clusters = parlay::filter(children, [&] (Cluster* x) { return !x->partner; });
                    if (adjusted_center_degree < 3) {
                        center->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                        if (adjusted_center_degree >= 1) {
                            Cluster* neighbor1 = center->get_neighbor();
                            if (!neighbor1->partner && neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                            if (adjusted_center_degree >= 2) {
                                Cluster* neighbor2 = center->get_other_neighbor(neighbor1);
                                if (!neighbor2->partner && neighbor2->parent == parent) local_root_clusters.push_back(center->get_other_neighbor(neighbor2));
                            }
                        }
                        center->insert_neighbors(children);
                        if (!center->partner) local_root_clusters.push_back(center);
                        parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                    }
                }
            }

            else if (max_degree == 2) {
                local_root_clusters.push_back(max);
                Cluster* neighbor1 = max->get_neighbor();
                Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                if (!neighbor1->partner && neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                else if (!neighbor2->partner && neighbor2->parent == parent) local_root_clusters.push_back(neighbor2);
                parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
            }

            else if (max_degree >= 3) {
                local_root_clusters = parlay::filter(children, [&] (Cluster* x) { return !x->partner; });
                if (max_degree - (children.size()-1) < 3) {
                    max->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                    Cluster* neighbor1 = max->get_neighbor();
                    if (!neighbor1->partner && neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                    if (max_degree - children.size() == 2) {
                        Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                        if (!neighbor2->partner && neighbor2->parent == parent) local_root_clusters.push_back(max->get_other_neighbor(neighbor2));
                    }
                    max->insert_neighbors(children);
                    parent->partner = (Cluster*) 1; // Use the partner field to mark a cluster for deletion
                } else {
                    auto iter = parlay::find(local_root_clusters, max);
                    if (iter != local_root_clusters.end()) local_root_clusters.erase(iter);
                }
            }

            parlay::parallel_for(0, children.size(), [&] (size_t i) {
                if (children[i]->partner) // should replace with parallel neighbor iteration
                    for (auto neighbor : children[i]->neighbors)
                        if (!neighbor->partner)
                            neighbor->delete_neighbor(children[i]);
            });

            parlay::parallel_for(0, local_root_clusters.size(), [&] (size_t i) {
                local_root_clusters[i]->parent = nullptr;
            });
            return local_root_clusters;
        }));

        parlay::parallel_for(0, del_clusters.size(), [&] (size_t i) {
            if (del_clusters[i]->partner) allocator::destroy(del_clusters[i]);
        });

        // ===============
        // PREP NEXT LEVEL
        // ===============
        root_clusters = parlay::append(next_root_clusters1, next_root_clusters2);
        del_clusters = std::move(next_del_clusters);
    }
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}
