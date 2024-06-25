#include "types.h"
#include "util.h"

/* This constant determines the maximum size of the vector of neighbors for each UFOCluster.
Any additional neighbors will be stored in the hash set for efficiency. */
#define UFO_ADJ_MAX 64

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    std::map<int, int> root_clusters_histogram;
#endif
// #define COLLECT_HEIGHT_STATS
#ifdef COLLECT_HEIGHT_STATS
    int max_height = 0;
#endif

long ufo_remove_ancestor_time = 0;
long ufo_recluster_tree_time = 0;


template<typename aug_t>
struct UFOCluster {
    // UFO cluster data
    std::vector<UFOCluster<aug_t>*> neighbors;
    std::unordered_set<UFOCluster<aug_t>*> neighbors_set;
    aug_t value; // Stores subtree values or cluster path values
    UFOCluster<aug_t>* parent = nullptr;
    // Constructor
    UFOCluster(aug_t value) : value(value) { neighbors.reserve(4); };
    // Helper functions
    int get_degree();
    bool high_degree();
    bool parent_high_fanout();
    bool contracts();
    bool contains_neighbor(UFOCluster<aug_t>* c);
    void insert_neighbor(UFOCluster<aug_t>* c, aug_t value = 1);
    void remove_neighbor(UFOCluster<aug_t>* c);
    UFOCluster<aug_t>* get_root();
};

template<typename aug_t>
class UFOTree {
public:
    // UFO tree interface
    UFOTree(vertex_t n, QueryType q = PATH, std::function<aug_t(aug_t, aug_t)> f = [](aug_t x, aug_t y)->aug_t{return x + y;}, aug_t id = 0, aug_t dval = 0);
    ~UFOTree();
    void link(vertex_t u, vertex_t v, aug_t value = 1);
    void cut(vertex_t u, vertex_t v);
    bool connected(vertex_t u, vertex_t v);
    // Testing helpers
    bool is_valid();
    int get_height(vertex_t v);
    void print_tree();
private:
    // Class data and parameters
    std::vector<UFOCluster<aug_t>> leaves;
    QueryType query_type;
    std::function<aug_t(aug_t, aug_t)> f;
    aug_t identity;
    aug_t default_value;
    std::vector<std::vector<UFOCluster<aug_t>*>> root_clusters;
    std::vector<std::pair<std::pair<UFOCluster<aug_t>*,UFOCluster<aug_t>*>,bool>> contractions;
    // lower_deg helps to identify clusters who became low degree during a deletion update
    std::vector<std::pair<UFOCluster<aug_t>*,int>> lower_deg;
    // Helper functions
    void remove_ancestors(UFOCluster<aug_t>* c, int start_level = 0);
    void recluster_tree();
    void disconnect_siblings(UFOCluster<aug_t>* c, int level);
    void insert_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v, aug_t value = 1);
    void remove_adjacency(UFOCluster<aug_t>* u, UFOCluster<aug_t>* v);
};

template<typename aug_t>
UFOTree<aug_t>::UFOTree(vertex_t n, QueryType q, std::function<aug_t(aug_t, aug_t)> f, aug_t id, aug_t d) :
query_type(q), f(f), identity(id), default_value(d) {
    leaves.resize(n, d);
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
}

template<typename aug_t>
UFOTree<aug_t>::~UFOTree() {
    #ifdef COLLECT_ROOT_CLUSTER_STATS
    std::cout << "Number of root clusters: Frequency" << std::endl;
        for (auto entry : root_clusters_histogram)
            std::cout << entry.first << "\t" << entry.second << std::endl;
    #endif
    #ifdef COLLECT_HEIGHT_STATS
        std::cout << "Maximum height of the tree: " << max_height << std::endl;
    #endif
    PRINT_TIMER("REMOVE ANCESTORS TIME", ufo_remove_ancestor_time);
    PRINT_TIMER("RECLUSTER TREE TIME", ufo_recluster_tree_time);
    return;
}

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
template<typename aug_t>
void UFOTree<aug_t>::link(vertex_t u, vertex_t v, aug_t value) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !connected(u,v));
    START_TIMER(ufo_remove_ancestor_timer);
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    STOP_TIMER(ufo_remove_ancestor_timer, ufo_remove_ancestor_time);
    insert_adjacency(&leaves[u], &leaves[v], value);
    START_TIMER(ufo_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(ufo_recluster_tree_timer, ufo_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
}

/* Cut vertex u and vertex v in the tree. */
template<typename aug_t>
void UFOTree<aug_t>::cut(vertex_t u, vertex_t v) {
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(leaves[u].contains_neighbor(&leaves[v]));
    auto curr_u = &leaves[u];
    auto curr_v = &leaves[v];
    while (curr_u != curr_v) {
        lower_deg.push_back({curr_u, curr_u->get_degree()-1});
        lower_deg.push_back({curr_v, curr_v->get_degree()-1});
        curr_u = curr_u->parent;
        curr_v = curr_v->parent;
    }
    START_TIMER(ufo_remove_ancestor_timer);
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    STOP_TIMER(ufo_remove_ancestor_timer, ufo_remove_ancestor_time);
    lower_deg.clear();
    remove_adjacency(&leaves[u], &leaves[v]);
    START_TIMER(ufo_recluster_tree_timer);
    recluster_tree();
    STOP_TIMER(ufo_recluster_tree_timer, ufo_recluster_tree_time);
    // Collect tree height stats at the end of each update
    #ifdef COLLECT_HEIGHT_STATS
        max_height = std::max(max_height, get_height(u));
        max_height = std::max(max_height, get_height(v));
    #endif
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
void UFOTree<aug_t>::remove_ancestors(UFOCluster<aug_t>* c, int start_level) {
    int level = start_level;
    auto prev = c;
    auto curr = c->parent;
    bool del = false;
    while (curr) {
        // Determine if curr is high degree or high fanout
        int curr_degree = curr->get_degree();
        int prev_degree = prev->get_degree();
        for (int i = 0; i < lower_deg.size(); i++) {
            if (lower_deg[i].first == curr) curr_degree = lower_deg[i].second;
            if (lower_deg[i].first == prev) prev_degree = lower_deg[i].second;
        }
        bool high_degree = (curr_degree > 2);
        bool high_fanout = false;
        if (prev->get_degree() == 1) {
            if (prev->neighbors[0]->parent == curr)
            if (prev->neighbors[0]->get_degree() - curr->get_degree() > 2) high_fanout = true;
        } else if (prev_degree - curr->get_degree() > 2) high_fanout = true;
        // Different cases for if curr will or will not be deleted later
        if (!high_degree && !high_fanout) { // We will delete curr next round
            disconnect_siblings(prev, level);
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry->remove_neighbor(prev); // Remove prev from adjacency
                auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
                if (position != root_clusters[level].end()) root_clusters[level].erase(position);
                delete prev;
            } else {
                prev->parent = nullptr;
                root_clusters[level].push_back(prev);
            }
            del = true;
        } else { // We will not delete curr next round
            if (del) { // Possibly delete prev
                for (auto entry : prev->neighbors)
                    entry->remove_neighbor(prev); // Remove prev from adjacency
                auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
                if (position != root_clusters[level].end()) root_clusters[level].erase(position);
                delete prev;
            } else if (prev->get_degree() <= 1) {
                prev->parent = nullptr;
                root_clusters[level].push_back(prev);
            }
            del = false;
        }
        // Update pointers
        prev = curr;
        curr = prev->parent;
        level++;
    }
    // DO LAST DELETIONS
    if (del) { // Possibly delete prev
        for (auto entry : prev->neighbors)
            entry->remove_neighbor(prev); // Remove prev from adjacency
        auto position = std::find(root_clusters[level].begin(), root_clusters[level].end(), prev);
        if (position != root_clusters[level].end()) root_clusters[level].erase(position);
        delete prev;
    } else root_clusters[level].push_back(prev);
}

template<typename aug_t>
void UFOTree<aug_t>::recluster_tree() {
    for (int level = 0; level < root_clusters.size(); level++) {
        if (root_clusters[level].empty()) continue;
        // Update root cluster stats if we are collecting them
        #ifdef COLLECT_ROOT_CLUSTER_STATS
            if (root_clusters_histogram.find(root_clusters[level].size()) == root_clusters_histogram.end())
                root_clusters_histogram[root_clusters[level].size()] = 1;
            else
                root_clusters_histogram[root_clusters[level].size()] += 1;
        #endif
        // Merge deg exactly 3 root clusters with all of its deg 1 neighbors
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 3) {
                auto parent = new UFOCluster<aug_t>(default_value);
                cluster->parent = parent;
                root_clusters[level+1].push_back(parent);
                bool first_contraction = true;
                for (auto entry : cluster->neighbors)
                    if (entry->get_degree() == 1) {
                        auto old_parent = entry->parent;
                        entry->parent = cluster->parent;
                        int lev = level+1;
                        while (old_parent) {
                            auto temp = old_parent;
                            old_parent = old_parent->parent;
                            auto position = std::find(root_clusters[lev].begin(), root_clusters[lev].end(), temp);
                            if (position != root_clusters[lev].end()) root_clusters[lev].erase(position);
                            delete temp;
                            lev++;
                        }
                        contractions.push_back({{entry, cluster},first_contraction});
                        first_contraction = false;
                    }
                if (first_contraction) contractions.push_back({{cluster, cluster}, true});
            }
        }
        // Always combine deg 1 root clusters with its neighboring cluster
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) {
                auto neighbor = cluster->neighbors[0];  // Deg 1 cluster only one neighbor
                if (neighbor->parent) {
                    if (neighbor->get_degree() == 2 && neighbor->contracts()) continue;
                    cluster->parent = neighbor->parent;
                    contractions.push_back({{cluster,neighbor},false});
                } else {
                    auto parent = new UFOCluster<aug_t>(default_value);
                    cluster->parent = parent;
                    neighbor->parent = parent;
                    root_clusters[level+1].push_back(parent);
                    contractions.push_back({{cluster,neighbor},true});
                }
                // For deg exactly 3 neighbor make sure all deg 1 neighbors combine with it
                if (neighbor->get_degree() == 3) {
                    for (auto entry : neighbor->neighbors) {
                        if (entry->get_degree() == 1 && entry->parent != neighbor->parent) {
                            auto old_parent = entry->parent;
                            entry->parent = neighbor->parent;
                            contractions.push_back({{entry, neighbor},false});
                            neighbor->parent->remove_neighbor(old_parent);
                            auto position = std::find(root_clusters[level+1].begin(), root_clusters[level+1].end(), old_parent);
                            if (position != root_clusters[level+1].end()) root_clusters[level+1].erase(position);
                            if (old_parent) delete old_parent;
                        }
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 2 root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto neighbor : cluster->neighbors) {
                    if (!neighbor->parent && (neighbor->get_degree() == 2)) {
                        auto parent = new UFOCluster<aug_t>(default_value);
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        root_clusters[level+1].push_back(parent);
                        contractions.push_back({{cluster,neighbor},true});
                        break;
                    }
                }
            }
        }
        // Combine deg 2 root clusters with deg 1 or 2 non-root clusters
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 2) {
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
                        if (neighbor->contracts()) continue;
                        cluster->parent = neighbor->parent;
                        contractions.push_back({{cluster,neighbor},false}); // The order here is important
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
                root_clusters[level+1].push_back(parent);
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
                    if (entry != c2) {
                        parent->insert_neighbor(entry->parent);
                        entry->parent->insert_neighbor(parent);
                    }
                }
                if (c1 != c2)
                for (auto entry : c2->neighbors) {
                    if (entry != c1) {
                        parent->insert_neighbor(entry->parent);
                        entry->parent->insert_neighbor(parent);
                    }
                }
            } else {
                // We ordered contractions so c1 would be the degree 1 cluster
                if (c1->get_degree() == 1) parent->remove_neighbor(c1->parent);
                // We ordered contractions so c2 is the one that had a parent already
                if (c1->get_degree() == 2) {
                    for (auto entry : c1->neighbors) if (entry != c2)
                        insert_adjacency(parent, entry->parent);
                }
                remove_ancestors(parent, level+1);
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
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

/* Helper function which takes a cluster c and the level of that cluster. The function
should find every cluster that shares a parent with c, disconnect it from their parent
and add it as a root cluster to be processed. */
template<typename aug_t>
void UFOTree<aug_t>::disconnect_siblings(UFOCluster<aug_t>* c, int level) {
    if (c->get_degree() == 1) {
        auto center = c->neighbors[0];
        if (center->parent == c->parent) {
            for (auto neighbor : center->neighbors) {
                if (neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            center->parent = nullptr;
            root_clusters[level].push_back(center);
        }
    } else {
        for (auto neighbor : c->neighbors) {
            if (neighbor->parent == c->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
    }
}

template<typename aug_t>
int UFOCluster<aug_t>::get_degree() {
    return neighbors.size() + neighbors_set.size();
}

/* Helper function which determines if the parent of a cluster is high fanout. This
cannot be determined from the parent itself since it does not store pointers to its
children, however the child stores adjacent clusters at its level. */
template<typename aug_t>
bool UFOCluster<aug_t>::parent_high_fanout() {
    assert(parent);
    if (get_degree() == 1) {
        if (neighbors[0]->parent == parent)
        if (neighbors[0]->get_degree() - parent->get_degree() > 2) return true;
    } else {
        if (get_degree() - parent->get_degree() > 2) return true;
    }
    return false;
}

// Helper function which returns whether this cluster combines with another cluster.
template<typename aug_t>
bool UFOCluster<aug_t>::contracts() {
    assert(get_degree() <= 3);
    for (auto neighbor : this->neighbors)
        if (neighbor->parent == this->parent)
            return true;
    return false;
}

template<typename aug_t>
bool UFOCluster<aug_t>::contains_neighbor(UFOCluster<aug_t>* c) {
    assert(get_degree() <= 3);
    if (std::find(neighbors.begin(), neighbors.end(), c) != neighbors.end()) return true;
    if (neighbors_set.find(c) != neighbors_set.end()) return true;
    return false;
}

template<typename aug_t>
void UFOCluster<aug_t>::insert_neighbor(UFOCluster<aug_t>* c, aug_t value) {
    assert(!contains_neighbor(c));
    if (contains_neighbor(c)) return;
    if (neighbors.size() < UFO_ADJ_MAX) neighbors.push_back(c);
    else neighbors_set.insert(c);
}

template<typename aug_t>
void UFOCluster<aug_t>::remove_neighbor(UFOCluster<aug_t>* c) {
    assert(contains_neighbor(c));
    auto position = std::find(neighbors.begin(), neighbors.end(), c);
    if (position != neighbors.end()) {
        if (neighbors_set.size() == 0) {
            std::iter_swap(position, neighbors.end()-1);
            neighbors.pop_back();
        } else { // Put an element from the set into the vector
            auto replacement = *neighbors_set.begin();
            *position = replacement;
            neighbors_set.erase(replacement);
        }
    } else neighbors_set.erase(c);
}

template<typename aug_t>
UFOCluster<aug_t>* UFOCluster<aug_t>::get_root() {
    UFOCluster<aug_t>* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}
