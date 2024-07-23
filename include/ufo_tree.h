#include "types.h"
#include "util.h"

/* These constants determines the maximum size of array of nieghbors and
the vector of neighbors for each UFOCluster. Any additional neighbors will
be stored in the hash set for efficiency. */
#define UFO_ARRAY_MAX 3
#define UFO_VECTOR_MAX 64

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
    UFOCluster<aug_t>* neighbors[UFO_ARRAY_MAX];
    std::vector<UFOCluster<aug_t>*> neighbors_vector;
    std::unordered_set<UFOCluster<aug_t>*> neighbors_set;
    aug_t value; // Stores subtree values or cluster path values
    UFOCluster<aug_t>* parent = nullptr;
    // Constructor
    UFOCluster(aug_t value) : neighbors(), parent(), value(value) {};
    // Helper functions
    int get_degree();
    bool high_degree();
    bool parent_high_fanout();
    bool contracts();
    bool contains_neighbor(UFOCluster<aug_t>* c);
    void insert_neighbor(UFOCluster<aug_t>* c, aug_t value = 1);
    void remove_neighbor(UFOCluster<aug_t>* c);
    UFOCluster<aug_t>* get_neighbor();
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
    // Clear all memory
    std::unordered_set<UFOCluster<aug_t>*> clusters;
    for (auto leaf : leaves) {
        auto curr = leaf.parent;
        while (curr) {
            clusters.insert(curr);
            curr = curr->parent;
        }
    }
    for (auto cluster : clusters) delete cluster;
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
            auto neighbor = prev->get_neighbor();
            if (neighbor->parent == curr)
            if (neighbor->get_degree() - curr->get_degree() > 2) high_fanout = true;
        } else if (prev_degree - curr_degree > 2) high_fanout = true;
        // Different cases for if curr will or will not be deleted later
        if (!high_degree && !high_fanout) { // We will delete curr next round
            disconnect_siblings(prev, level);
            if (del) { // Possibly delete prev
                for (auto neighbor : prev->neighbors)
                    if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
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
                for (auto neighbor : prev->neighbors)
                    if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
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
        for (auto neighbor : prev->neighbors)
            if (neighbor) neighbor->remove_neighbor(prev); // Remove prev from adjacency
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
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto old_parent = neighbor->parent;
                        neighbor->parent = cluster->parent;
                        int lev = level+1;
                        while (old_parent) {
                            auto temp = old_parent;
                            old_parent = old_parent->parent;
                            auto position = std::find(root_clusters[lev].begin(), root_clusters[lev].end(), temp);
                            if (position != root_clusters[lev].end()) root_clusters[lev].erase(position);
                            delete temp;
                            lev++;
                        }
                        if (first_contraction) contractions.push_back({{neighbor, cluster},true});
                        first_contraction = false;
                    }
                }
                if (first_contraction) contractions.push_back({{cluster, cluster}, true});
            }
        }
        // Always combine deg 1 root clusters with its neighboring cluster
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() == 1) {
                auto neighbor = cluster->get_neighbor();  // Deg 1 cluster only one neighbor
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
                        if (entry && entry->get_degree() == 1 && entry->parent != neighbor->parent) {
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
                    if (neighbor && !neighbor->parent && (neighbor->get_degree() == 2)) {
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
                    if (neighbor && neighbor->parent && (neighbor->get_degree() == 1 || neighbor->get_degree() == 2)) {
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
                for (int i = 0; i < 3; ++i) parent->neighbors[i] = nullptr;
                for (auto neighbor : c1->neighbors) {
                    if (neighbor && neighbor != c2 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
                for (auto neighbor : c1->neighbors_vector) {
                    if (neighbor != c2 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
                if (c1 != c2)
                for (auto neighbor : c2->neighbors) {
                    if (neighbor && neighbor != c1 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
                for (auto neighbor : c2->neighbors_vector) {
                    if (neighbor != c1 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
            } else {
                // We ordered contractions so c2 is the one that had a parent already
                if (c1->get_degree() == 2) {
                    for (auto neighbor : c1->neighbors) if (neighbor && neighbor != c2)
                        insert_adjacency(parent, neighbor->parent);
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
        auto center = c->get_neighbor();
        if (c->parent != center->parent) return;
        assert(center->get_degree() <= 5);
        if (center->parent == c->parent) {
            for (auto neighbor : center->neighbors) {
                if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            for (auto neighbor : center->neighbors_vector) {
                if (neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            center->parent = nullptr;
            root_clusters[level].push_back(center);
        }
    } else {
        assert(c->get_degree() <= 5);
        for (auto neighbor : c->neighbors) {
            if (neighbor && neighbor->parent == c->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
        for (auto neighbor : c->neighbors_vector) {
            if (neighbor->parent == c->parent) {
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
    }
}

template<typename aug_t>
int UFOCluster<aug_t>::get_degree() {
    int deg = 0;
    for (auto neighbor : this->neighbors) if (neighbor) deg++;
    return deg + neighbors_vector.size() + neighbors_set.size();
}

/* Helper function which determines if the parent of a cluster is high fanout. This
cannot be determined from the parent itself since it does not store pointers to its
children, however the child stores adjacent clusters at its level. */
template<typename aug_t>
bool UFOCluster<aug_t>::parent_high_fanout() {
    assert(parent);
    if (get_degree() == 1) {
        auto neighbor = get_neighbor();
        if (neighbor->parent == parent)
        if (neighbor->get_degree() - parent->get_degree() > 2) return true;
    } else {
        if (get_degree() - parent->get_degree() > 2) return true;
    }
    return false;
}

// Helper function which returns whether this cluster combines with another cluster.
template<typename aug_t>
bool UFOCluster<aug_t>::contracts() {
    assert(get_degree() <= 3);
    for (auto neighbor : this->neighbors) if (neighbor && neighbor->parent == this->parent) return true;
    return false;
}

template<typename aug_t>
bool UFOCluster<aug_t>::contains_neighbor(UFOCluster<aug_t>* c) {
    // This should only be used in testing
    for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
    for (auto neighbor : neighbors_vector) if (neighbor == c) return true;
    if (neighbors_set.find(c) != neighbors_set.end()) return true;
    return false;
}

template<typename aug_t>
void UFOCluster<aug_t>::insert_neighbor(UFOCluster<aug_t>* c, aug_t value) {
    for (int i = 0; i < 3; ++i) if (this->neighbors[i] == c) return;
    assert(!contains_neighbor(c));
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == nullptr) {
            this->neighbors[i] = c;
            return;
        }
    }
    if (neighbors_vector.size() < UFO_VECTOR_MAX) neighbors_vector.push_back(c);
    else neighbors_set.insert(c);
}

template<typename aug_t>
void UFOCluster<aug_t>::remove_neighbor(UFOCluster<aug_t>* c) {
    assert(contains_neighbor(c));
    for (int i = 0; i < 3; ++i) {
        if (this->neighbors[i] == c) {
            this->neighbors[i] = nullptr;
            if (!neighbors_set.empty()) { // Put an element from the set into the array
                auto replacement = *neighbors_set.begin();
                this->neighbors[i] = replacement;
                neighbors_set.erase(replacement);
            } else if (!neighbors_vector.empty()) { // Put an element from the vector into the array
                auto replacement = neighbors_vector.back();
                this->neighbors[i] = replacement;
                neighbors_vector.pop_back();
            }
            return;
        }
    }
    auto position = std::find(neighbors_vector.begin(), neighbors_vector.end(), c);
    if (position != neighbors_vector.end()) {
        if (neighbors_set.empty()) {
            std::iter_swap(position, neighbors_vector.end()-1);
            neighbors_vector.pop_back();
        } else { // Put an element from the set into the vector
            auto replacement = *neighbors_set.begin();
            *position = replacement;
            neighbors_set.erase(replacement);
        }
    } else neighbors_set.erase(c);
}

template<typename aug_t>
UFOCluster<aug_t>* UFOCluster<aug_t>::get_neighbor() {
    for (auto neighbor : neighbors)
        if (neighbor) return neighbor;
    std::cerr << "No neighbor found to return." << std::endl;
    std::abort();
}

template<typename aug_t>
UFOCluster<aug_t>* UFOCluster<aug_t>::get_root() {
    UFOCluster<aug_t>* curr = this;
    while (curr->parent) curr = curr->parent;
    return curr;
}
