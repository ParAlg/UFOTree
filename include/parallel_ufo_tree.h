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
    void recluster_tree(parlay::sequence<std::pair<Cluster*, Cluster*>>& updates, UpdateType update_type);
};

template <typename aug_t>
ParallelUFOTree<aug_t>::ParallelUFOTree(vertex_t n, vertex_t k) : leaves(n) {}

template <typename aug_t>
ParallelUFOTree<aug_t>::~ParallelUFOTree() {
    allocator::finish();
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_link(parlay::sequence<std::pair<int, int>>& links) {
    auto updates = parlay::tabulate(2*links.size(), [&] (size_t i) {
        if (i % 2 == 0) return std::make_pair(&leaves[links[i].first], &leaves[links[i].second]);
        return std::make_pair(&leaves[links[i].second], &leaves[links[i].first]);
    });
    recluster_tree(updates, INSERT);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::batch_cut(parlay::sequence<std::pair<int, int>>& cuts) {
    auto updates = parlay::tabulate(2*cuts.size(), [&] (size_t i) {
        if (i % 2 == 0) return std::make_pair(&leaves[cuts[i].first], &leaves[cuts[i].second]);
        return std::make_pair(&leaves[cuts[i].second], &leaves[cuts[i].first]);
    });
    recluster_tree(updates, DELETE);
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::recluster_tree(parlay::sequence<std::pair<ParallelUFOCluster<aug_t>*, ParallelUFOCluster<aug_t>*>>& updates, UpdateType update_type) {
    parlay::sequence<Cluster*> root_clusters;
    parlay::sequence<Cluster*> del_clusters;
    
    auto endpoint_groups = parlay::group_by_key(updates);
    auto parents = parlay::delayed_tabulate(endpoint_groups.size(), [&] (size_t i) {
        return std::make_pair(endpoint_groups[i].first->parent, endpoint_groups[i].first);
    });
    auto parent_groups = parlay::group_by_key(parents);

    // The intial del clusters are just all parents of vertices that got an update.
    del_clusters = parlay::map(parent_groups, [&] (auto group) {
        return group.first;
    });

    /* Determine the root clusters by figuring out which parents will actually get deleted. Any children of
    deleted level i+1 clusters should be root clusters. Additionally, any low degree cluster will become a
    root cluster regardless of whether its parent is deleted. Also unset the parent pointer of root clusters. */
    root_clusters = parlay::flatten(parlay::delayed_tabulate(parent_groups.size()/* +1 */, [&] (size_t i) {
        // if (i == parent_groups.size()) return next_root_clusters;
        auto& [parent, children] = parent_groups[i];
        parlay::sequence<Cluster*> local_root_clusters;
        if (parent == nullptr) return local_root_clusters;

        Cluster* max = *parlay::max_element(children, [&] (Cluster* x, Cluster* y) { return x->get_degree() < y->get_degree(); });
        size_t max_degree = max->get_degree();

        if (max_degree == 1) {
            local_root_clusters = children; // We can probably pass by reference here
            Cluster* center = max->get_neighbor();
            size_t center_degree = center->get_degree();
            if (center_degree - children.size() < 3) {
                center->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                Cluster* neighbor1 = center->get_neighbor();
                if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                if (center_degree - children.size() == 2) {
                    Cluster* neighbor2 = center->get_other_neighbor(neighbor1);
                    if (neighbor2->parent == parent) local_root_clusters.push_back(center->get_other_neighbor(neighbor2));
                }
                center->insert_neighbors(children);
            }
        }

        else if (max_degree == 2) {
            local_root_clusters.push_back(max);
            Cluster* neighbor = max->get_neighbor();
            if (neighbor->parent = parent) local_root_clusters.push_back(neighbor);
            else local_root_clusters.push_back(max->get_other_neighbor(neighbor));
        }

        else if (max_degree >= 3) {
            local_root_clusters = children; // We can probably pass by reference here
            if (max_degree - (children.size()-1) < 3) {
                max->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                Cluster* neighbor1 = max->get_neighbor();
                if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                if (max_degree - children.size() == 2) {
                    Cluster* neighbor2 = max->get_other_neighbor(neighbor1);
                    if (neighbor2->parent == parent) local_root_clusters.push_back(max->get_other_neighbor(neighbor2));
                }
                max->insert_neighbors(children);
            } else {
                local_root_clusters.erase(parlay::find(local_root_clusters, max));
            }
        }

        return local_root_clusters;
    }));

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
        // PHASE 1
        // =======
        auto parents = parlay::delayed_tabulate(del_clusters.size(), [&] (size_t i) {
            return std::make_pair(del_clusters[i]->parent, del_clusters[i]);
        });
        auto parent_groups = parlay::group_by_key(parents);
        // Get the next del clusters
        auto next_del_clusters = parlay::delayed::map_maybe(parent_groups, [&] (auto group) -> std::optional<Cluster*> {
            if (!group.first) return std::nullopt;
            return group.first;
        });
        // Get next root clusters, clear their parent pointers, delete some current del clusters
        auto next_root_clusters = parlay::delayed::flatten(parlay::delayed_tabulate(parent_groups.size(), [&] (size_t i) {
            auto& [parent, children] = parent_groups[i];
            parlay::sequence<Cluster*> local_root_clusters;
            if (parent == nullptr) return local_root_clusters;

            Cluster* max = *parlay::max_element(children, [&] (Cluster* x, Cluster* y) { return x->get_degree() < y->get_degree(); });
            size_t max_degree = max->get_degree();

            if (max_degree == 1) {
                local_root_clusters = std::move(children);
                Cluster* center = max->get_neighbor();
                size_t center_degree = center->get_degree();
                if (center_degree - children.size() < 3) {
                    center->delete_neighbors(children); // Should replace delete / reinsert with parallel iteration
                    Cluster* neighbor1 = center->get_neighbor();
                    if (neighbor1->parent == parent) local_root_clusters.push_back(neighbor1);
                    if (center_degree - children.size() == 2) {
                        Cluster* neighbor2 = center->get_other_neighbor(neighbor1);
                        if (neighbor2->parent == parent) local_root_clusters.push_back(center->get_other_neighbor(neighbor2));
                    }
                    center->insert_neighbors(children);
                    allocator::destroy(parent);
                }
            }

            else if (max_degree == 2) {
                local_root_clusters.push_back(max);
                Cluster* neighbor = max->get_neighbor();
                if (neighbor->parent = parent) local_root_clusters.push_back(neighbor);
                else local_root_clusters.push_back(max->get_other_neighbor(neighbor));
                allocator::destroy(parent);
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
                    allocator::destroy(parent);
                } else {
                    local_root_clusters.erase(parlay::find(local_root_clusters, max));
                }
            }

            parlay::parallel_for(0, local_root_clusters.size(), [&] (size_t i) {
                local_root_clusters[i]->parent = nullptr;
            });
            return local_root_clusters;
        }));
        // =======
        // PHASE 2
        // =======

        // insert/delete the edges at this level
        parlay::parallel_for(0, endpoint_groups.size(), [&] (size_t i) {
            auto& [vertex, neighbors] = endpoint_groups[i];
            if (update_type == INSERT) vertex->insert_neighbors(neighbors);
            else vertex->delete_neighbors(neighbors);
        });

        // Prepare next level
        updates = parlay::filter(updates, [&] (std::pair<Cluster*, Cluster*> x) {
            return x.first->parent != x.second->parent;
        });

        /* We also need to handle the case where a root cluster became high degree because of many
        insertions, but it needs to find a degree 1 neighbor that is a non-root cluster because it
        previously happened to not contract with this cluster. Fortunately we can detect when this
        will happen by looking at the number of edges inserted, and we can grab that deg 1 neighbor
        easily before inserting all the edges and making this cluster high degree.
        This might mean populating the adjacency lists should happen during the phase 3 of the next
        level up. */
    }
}

template <typename aug_t>
bool ParallelUFOTree<aug_t>::connected(vertex_t u, vertex_t v) {
    return (leaves[u].get_root() == leaves[v].get_root());
}

}
