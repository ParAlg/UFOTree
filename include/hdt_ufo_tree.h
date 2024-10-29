#pragma once
#include "types.h"
#include "util.h"
#include <absl/container/flat_hash_set.h>
#include <absl/container/flat_hash_map.h>

/* These constants determines the maximum size of array of nieghbors and
the vector of neighbors for each HDTUFOCluster. Any additional neighbors will
be stored in the hash set for efficiency. */
#define HDTUFO_ARRAY_MAX 3

// #define COLLECT_ROOT_CLUSTER_STATS
#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> root_clusters_histogram;
#endif

struct HDTUFOCluster {
    // HDTUFO cluster data
    HDTUFOCluster* parent = nullptr;
    HDTUFOCluster* neighbors[HDTUFO_ARRAY_MAX];
    absl::flat_hash_set<HDTUFOCluster*>* neighbors_set = nullptr;
    vertex_t size = 1;
    vertex_t vertex_mark = NONE;   // store the id of a marked vertex in the cluster
    vertex_t edge_mark = NONE;     // store the id of a vertex in the cluster with a marked incident edge
    // Keeps track of siblings that are marked and share a parent with this cluster
    absl::flat_hash_set<HDTUFOCluster*>* vertex_marked_siblings;
    absl::flat_hash_set<HDTUFOCluster*>* edge_marked_siblings;
    // Constructor
    HDTUFOCluster() :  parent(), neighbors(), size(1), vertex_mark(NONE), edge_mark(NONE), vertex_marked_siblings(), edge_marked_siblings(){};
    // Helper functions
    bool contracts() {
        assert(get_degree() <= 3);
        for (auto neighbor : this->neighbors) if (neighbor && neighbor->parent == this->parent) return true;
        return false;
    }
    HDTUFOCluster* get_neighbor() {
        for (auto neighbor : neighbors)
            if (neighbor) return neighbor;
        std::cerr << "No neighbor found to return." << std::endl;
        std::abort();
    }
    vertex_t get_sibling_vertex_mark() {
        // for (auto neighbor : neighbors)
        //     if (neighbor && neighbor->parent == parent && neighbor->vertex_mark != NONE) return neighbor->vertex_mark;
        if (!vertex_marked_siblings->empty()) return (*vertex_marked_siblings->begin())->vertex_mark;
        return this->vertex_mark;
    }
    vertex_t get_sibling_edge_mark() {
        // for (auto neighbor : neighbors)
        //     if (neighbor && neighbor->parent == parent && neighbor->edge_mark != NONE) return neighbor->edge_mark;
        if (!edge_marked_siblings->empty()) return (*edge_marked_siblings->begin())->edge_mark;
        return this->edge_mark;
    }
    HDTUFOCluster* get_root() {
        HDTUFOCluster* curr = this;
        while (curr->parent) curr = curr->parent;
        return curr;
    }
    int get_degree() {
        int deg = 0;
        for (auto neighbor : this->neighbors) if (neighbor) deg++;
        if (neighbors_set) deg += neighbors_set->size();
        return deg;
    }
    bool parent_high_fanout() {
        assert(parent);
        int parent_degree = parent->get_degree();
        if (get_degree() == 1) {
            auto neighbor = get_neighbor();
            if (neighbor->parent == parent)
            if (neighbor->get_degree() - parent_degree > 2) return true;
        } else {
            if (get_degree() - parent_degree > 2) return true;
        }
        return false;
    }
    bool contains_neighbor(HDTUFOCluster* c) {
        for (auto neighbor : neighbors) if (neighbor && neighbor == c) return true;
        if (neighbors_set && neighbors_set->find(c) != neighbors_set->end()) return true;
        return false;
    }
    void insert_neighbor(HDTUFOCluster* c) {
        for (int i = 0; i < HDTUFO_ARRAY_MAX; ++i) if (this->neighbors[i] == c) return;
        assert(!contains_neighbor(c));
        for (int i = 0; i < HDTUFO_ARRAY_MAX; ++i) {
            if (this->neighbors[i] == nullptr) {
                this->neighbors[i] = c;
                return;
            }
        }
        if (!neighbors_set)
            neighbors_set = new absl::flat_hash_set<HDTUFOCluster*>();
        neighbors_set->insert(c);
    }
    void remove_neighbor(HDTUFOCluster* c) {
        assert(contains_neighbor(c));
        for (int i = 0; i < 3; ++i) {
            if (this->neighbors[i] == c) {
                this->neighbors[i] = nullptr;
                if (neighbors_set) { // Put an element from the set into the array
                    auto replacement = *neighbors_set->begin();
                    this->neighbors[i] = replacement;
                    neighbors_set->erase(replacement);
                    vertex_marked_siblings->erase(replacement);
                    edge_marked_siblings->erase(replacement);
                    if (neighbors_set->empty()) {
                        delete neighbors_set;
                        neighbors_set = nullptr;
                    }
                    if (vertex_marked_siblings->empty()) {
                        delete vertex_marked_siblings;
                        vertex_marked_siblings = nullptr;
                    }
                    if (edge_marked_siblings->empty()) {
                        delete edge_marked_siblings;
                        edge_marked_siblings = nullptr;
                    }
                }
                return;
            }
        }
        neighbors_set->erase(c);
        vertex_marked_siblings->erase(c);
        edge_marked_siblings->erase(c);
        if (neighbors_set->empty()) {
            delete neighbors_set;
            neighbors_set = nullptr;
        }
        if (vertex_marked_siblings->empty()) {
            delete vertex_marked_siblings;
            vertex_marked_siblings = nullptr;
        }
        if (edge_marked_siblings->empty()) {
            delete edge_marked_siblings;
            edge_marked_siblings = nullptr;
        }
    }
    size_t calculate_size(){
        size_t memory = sizeof(HDTUFOCluster);
        if (neighbors_set)
            memory += neighbors_set->bucket_count() * sizeof(HDTUFOCluster*);
        if (vertex_marked_siblings)
            memory += vertex_marked_siblings->bucket_count() * sizeof(HDTUFOCluster*);
        if (edge_marked_siblings)
            memory += edge_marked_siblings->bucket_count() * sizeof(HDTUFOCluster*);
        return memory;
    }
};

class HDTUFOTree {
using Cluster = HDTUFOCluster;
public:
    // HDT Dynamic Tree interface
    HDTUFOTree(vertex_t n);
    ~HDTUFOTree();
    void AddEdge(edge_t e);
    void DeleteEdge(edge_t e);
    void MarkVertex(vertex_t v, bool mark); // Marked vertices have an incident level i non-tree edge
    void MarkEdge(edge_t e, bool mark);     // Marked edges represent level i tree edges (not all tree edges)
    std::optional<vertex_t> GetMarkedVertexInTree(vertex_t v);
    std::optional<edge_t> GetMarkedEdgeInTree(vertex_t v);
    vertex_t GetSizeOfTree(vertex_t v);
    bool IsConnected(vertex_t u, vertex_t v);
    // Testing helpers
    size_t space();
    size_t count_nodes();
    size_t get_height();
    bool is_valid();
    void print_tree();
private:
    // Class data and parameters
    std::vector<Cluster> leaves;
    std::vector<absl::flat_hash_set<vertex_t>*> marked_tree_edges;
    std::vector<std::vector<Cluster*>> root_clusters;
    std::vector<std::pair<std::pair<Cluster*,Cluster*>,bool>> contractions;
    std::vector<std::pair<Cluster*,int>> lower_deg; // lower_deg helps to identify clusters who became low degree during a deletion update
    // Helper functions
    void remove_ancestors(Cluster* c, int start_level = 0);
    void recluster_tree();
    void disconnect_siblings(Cluster* c, int level);
    void insert_adjacency(Cluster* u, Cluster* v);
    void remove_adjacency(Cluster* u, Cluster* v);
    void add_vertex_mark(vertex_t v);
    void add_edge_mark(vertex_t v);
    void remove_vertex_mark(vertex_t v);
    void remove_edge_mark(vertex_t v);
    void recompute_parent_vertex_mark(Cluster* c);
    void recompute_parent_edge_mark(Cluster* c);
    void add_sibling_marks(Cluster* c1, Cluster* c2);
};

HDTUFOTree::HDTUFOTree(vertex_t n) {
    leaves.resize(n);
    marked_tree_edges.resize(n);
    for (int i = 0; i < n; ++i) marked_tree_edges[i] = new absl::flat_hash_set<vertex_t>;
    root_clusters.resize(max_tree_height(n));
    contractions.reserve(12);
}

HDTUFOTree::~HDTUFOTree() {
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
    #ifdef COLLECT_ROOT_CLUSTER_STATS
    std::cout << "Number of root clusters: Frequency" << std::endl;
        for (auto entry : root_clusters_histogram)
            std::cout << entry.first << "\t" << entry.second << std::endl;
    #endif
}

size_t HDTUFOTree::space() {
    std::unordered_set<Cluster*> visited;
    size_t memory = sizeof(HDTUFOTree);
    for(auto cluster : leaves){
        memory += cluster.calculate_size();
        auto parent = cluster.parent;
        while(parent != nullptr && visited.count(parent) == 0){
            memory += parent->calculate_size();
            visited.insert(parent);
            parent = parent->parent;
        }
    }
    return memory;
}

size_t HDTUFOTree::count_nodes() {
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

size_t HDTUFOTree::get_height() {
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

/* Link vertex u and vertex v in the tree. Optionally include an
augmented value for the new edge (u,v). If no augmented value is
provided, the default value is 1. */
void HDTUFOTree::AddEdge(edge_t e) {
    vertex_t u = e.first;
    vertex_t v = e.second;
    assert(u >= 0 && u < leaves.size() && v >= 0 && v < leaves.size());
    assert(u != v && !IsConnected(u,v));
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    insert_adjacency(&leaves[u], &leaves[v]);
    recluster_tree();
}

/* Cut vertex u and vertex v in the tree. */
void HDTUFOTree::DeleteEdge(edge_t e) {
    vertex_t u = e.first;
    vertex_t v = e.second;
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
    remove_ancestors(&leaves[u]);
    remove_ancestors(&leaves[v]);
    lower_deg.clear();
    remove_adjacency(&leaves[u], &leaves[v]);
    recluster_tree();
    MarkEdge(e, false);
}

/* Removes the ancestors of cluster c that are not high degree nor
high fan-out and add them to root_clusters. */
void HDTUFOTree::remove_ancestors(Cluster* c, int start_level) {
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
                prev->parent->size -= prev->size;
                recompute_parent_vertex_mark(prev);
                recompute_parent_edge_mark(prev);
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
                prev->parent->size -= prev->size;
                recompute_parent_vertex_mark(prev);
                recompute_parent_edge_mark(prev);
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

void HDTUFOTree::recluster_tree() {
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
                auto parent = new Cluster();
                cluster->parent = parent;
                parent->size = cluster->size;
                parent->vertex_mark = cluster->vertex_mark;
                parent->edge_mark = cluster->edge_mark;
                root_clusters[level+1].push_back(parent);
                bool first_contraction = true;
                for (auto neighbor : cluster->neighbors) {
                    if (neighbor && neighbor->get_degree() == 1) {
                        auto old_parent = neighbor->parent;
                        neighbor->parent = cluster->parent;
                        cluster->parent->size += neighbor->size;
                        if (cluster->vertex_mark != NONE) parent->vertex_mark = cluster->vertex_mark;
                        if (cluster->edge_mark != NONE) parent->edge_mark = cluster->edge_mark;
                        add_sibling_marks(cluster, neighbor);
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
                    neighbor->parent->size += cluster->size;
                    if (cluster->vertex_mark != NONE) neighbor->parent->vertex_mark = cluster->vertex_mark;
                    if (cluster->edge_mark != NONE) neighbor->parent->edge_mark = cluster->edge_mark;
                    add_sibling_marks(cluster, neighbor);
                    contractions.push_back({{cluster,neighbor},false});
                } else {
                    auto parent = new Cluster();
                    cluster->parent = parent;
                    neighbor->parent = parent;
                    parent->size = cluster->size + neighbor->size;
                    if (cluster->vertex_mark != NONE) parent->vertex_mark = cluster->vertex_mark;
                    if (cluster->edge_mark != NONE) parent->edge_mark = cluster->edge_mark;
                    if (neighbor->vertex_mark != NONE) parent->vertex_mark = neighbor->vertex_mark;
                    if (neighbor->edge_mark != NONE) parent->edge_mark = neighbor->edge_mark;
                    add_sibling_marks(cluster, neighbor);
                    root_clusters[level+1].push_back(parent);
                    contractions.push_back({{cluster,neighbor},true});
                }
                // For deg exactly 3 neighbor make sure all deg 1 neighbors combine with it
                if (neighbor->get_degree() == 3) {
                    for (auto entry : neighbor->neighbors) {
                        if (entry && entry->get_degree() == 1 && entry->parent != neighbor->parent) {
                            auto old_parent = entry->parent;
                            entry->parent = neighbor->parent;
                            neighbor->parent->size += entry->size;
                            if (entry->vertex_mark != NONE) neighbor->parent->vertex_mark = entry->vertex_mark;
                            if (entry->edge_mark != NONE) neighbor->parent->edge_mark = entry->edge_mark;
                            add_sibling_marks(neighbor, entry);
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
                        auto parent = new Cluster();
                        cluster->parent = parent;
                        neighbor->parent = parent;
                        parent->size = cluster->size + neighbor->size;
                        if (cluster->vertex_mark != NONE) parent->vertex_mark = cluster->vertex_mark;
                        if (cluster->edge_mark != NONE) parent->edge_mark = cluster->edge_mark;
                        if (neighbor->vertex_mark != NONE) parent->vertex_mark = neighbor->vertex_mark;
                        if (neighbor->edge_mark != NONE) parent->edge_mark = neighbor->edge_mark;
                        add_sibling_marks(cluster, neighbor);
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
                        neighbor->parent->size += cluster->size;
                        if (cluster->vertex_mark != NONE) neighbor->parent->vertex_mark = cluster->vertex_mark;
                        if (cluster->edge_mark != NONE) neighbor->parent->edge_mark = cluster->edge_mark;
                        add_sibling_marks(cluster, neighbor);
                        contractions.push_back({{cluster,neighbor},false}); // The order here is important
                        break;
                    }
                }
            }
        }
        // Add remaining uncombined clusters to the next level
        for (auto cluster : root_clusters[level]) {
            if (!cluster->parent && cluster->get_degree() > 0) {
                auto parent = new Cluster();
                cluster->parent = parent;
                parent->size = cluster->size;
                parent->vertex_mark = cluster->vertex_mark;
                parent->edge_mark = cluster->edge_mark;
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
                for (int i = 0; i < HDTUFO_ARRAY_MAX; i++) {
                    auto neighbor = c1->neighbors[i];
                    if (neighbor && neighbor != c2 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
                if (c1->neighbors_set)
                for (auto neighbor : *c1->neighbors_set) {
                    if (neighbor && neighbor != c2 && parent != neighbor->parent) {
                        parent->insert_neighbor(neighbor->parent);
                        neighbor->parent->insert_neighbor(parent);
                    }
                }
                if (c1 != c2) {
                    for (int i = 0; i < HDTUFO_ARRAY_MAX; i++) {
                        auto neighbor = c2->neighbors[i];
                        if (neighbor && neighbor != c1 && parent != neighbor->parent) {
                            parent->insert_neighbor(neighbor->parent);
                            neighbor->parent->insert_neighbor(parent);
                        }
                    }
                    if  (c2->neighbors_set)
                    for (auto neighbor : *c2->neighbors_set){
                        if (neighbor && neighbor != c1 && parent != neighbor->parent) {
                            parent->insert_neighbor(neighbor->parent);
                            neighbor->parent->insert_neighbor(parent);
                        }
                    }
                }
            } else {
                // We ordered contractions so c2 is the one that had a parent already
                if (c1->get_degree() == 2) {
                    for (int i = 0; i < HDTUFO_ARRAY_MAX; i++){
                        auto neighbor = c1->neighbors[i];
                        if (neighbor && neighbor != c2)
                        insert_adjacency(parent, neighbor->parent);
                    }
                }
                remove_ancestors(parent, level+1);
            }
        }
        // Clear the contents of this level
        root_clusters[level].clear();
        contractions.clear();
    }
}

void HDTUFOTree::insert_adjacency(Cluster* u, Cluster* v) {
    auto curr_u = u;
    auto curr_v = v;
    while (curr_u && curr_v && curr_u != curr_v) {
        curr_u->insert_neighbor(curr_v);
        curr_v->insert_neighbor(curr_u);
        curr_u = curr_u->parent;
        curr_v = curr_v->parent;
    }
}

void HDTUFOTree::remove_adjacency(Cluster* u, Cluster* v) {
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
void HDTUFOTree::disconnect_siblings(Cluster* c, int level) {
    if (c->get_degree() == 1) {
        auto center = c->get_neighbor();
        if (c->parent != center->parent) return;
        assert(center->get_degree() <= 5);
        if (center->parent == c->parent) {
            for (auto neighbor : center->neighbors) {
                if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent->size -= neighbor->size;
                    recompute_parent_vertex_mark(neighbor);
                    recompute_parent_edge_mark(neighbor);
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            if (center->neighbors_set)
            for (auto neighbor : *center->neighbors_set) {
                if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                    neighbor->parent->size -= neighbor->size;
                    recompute_parent_vertex_mark(neighbor);
                    recompute_parent_edge_mark(neighbor);
                    neighbor->parent = nullptr; // Set sibling parent pointer to null
                    root_clusters[level].push_back(neighbor); // Keep track of root clusters
                }
            }
            center->parent->size -= center->size;
            recompute_parent_vertex_mark(center);
            recompute_parent_edge_mark(center);
            center->parent = nullptr;
            root_clusters[level].push_back(center);
        }
    } else {
        assert(c->get_degree() <= 5);
        for (auto neighbor : c->neighbors) {
            if (neighbor && neighbor->parent == c->parent) {
                neighbor->parent->size -= neighbor->size;
                recompute_parent_vertex_mark(neighbor);
                recompute_parent_edge_mark(neighbor);
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
        if (c->neighbors_set)
        for (auto neighbor : *c->neighbors_set) {
            if (neighbor && neighbor->parent == c->parent && neighbor != c) {
                neighbor->parent->size -= neighbor->size;
                recompute_parent_vertex_mark(neighbor);
                recompute_parent_edge_mark(neighbor);
                neighbor->parent = nullptr; // Set sibling parent pointer to null
                root_clusters[level].push_back(neighbor); // Keep track of root clusters
            }
        }
    }
}

void HDTUFOTree::MarkVertex(vertex_t v, bool mark) {
    if (mark) add_vertex_mark(v);
    else remove_vertex_mark(v);
}

void HDTUFOTree::MarkEdge(edge_t e, bool mark) {
    if (mark) {
        if (marked_tree_edges[e.first]->empty()) add_edge_mark(e.first);
        marked_tree_edges[e.first]->insert(e.second);
        if (marked_tree_edges[e.second]->empty()) add_edge_mark(e.second);
        marked_tree_edges[e.second]->insert(e.first);
    }
    else {
        marked_tree_edges[e.first]->erase(e.second);
        if (marked_tree_edges[e.first]->empty()) remove_edge_mark(e.first);
        marked_tree_edges[e.second]->erase(e.first);
        if (marked_tree_edges[e.second]->empty()) remove_edge_mark(e.second);
    }
}

void HDTUFOTree::add_vertex_mark(vertex_t v) {
    auto curr = &leaves[v];
    while (curr) {
        if (curr->vertex_mark == NONE) {
            for (auto neighbor : curr->neighbors) {
                if (neighbor && neighbor->neighbors_set && neighbor->neighbors_set->contains(curr)) {
                    if (!neighbor->vertex_marked_siblings)
                        neighbor->vertex_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
                    neighbor->vertex_marked_siblings->insert(curr);
                }
            }
        }
        curr->vertex_mark = v;
        curr = curr->parent;
    }
}

void HDTUFOTree::add_edge_mark(vertex_t v) {
    auto curr = &leaves[v];
    while (curr) {
        if (curr->edge_mark == NONE) {
            for (auto neighbor : curr->neighbors) {
                if (neighbor && neighbor->neighbors_set && neighbor->neighbors_set->contains(curr)) {
                    if (!neighbor->edge_marked_siblings)
                        neighbor->edge_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
                    neighbor->edge_marked_siblings->insert(curr);
                }
            }
        }
        curr->edge_mark = v;
        curr = curr->parent;
    }
}

void HDTUFOTree::remove_vertex_mark(vertex_t v) {
    auto curr = &leaves[v];
    curr->vertex_mark = NONE;
    for (auto neighbor : curr->neighbors) {
        if (neighbor && neighbor->neighbors_set && neighbor->neighbors_set->contains(curr)) {
            if (!neighbor->vertex_marked_siblings)
                neighbor->vertex_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
            neighbor->vertex_marked_siblings->insert(curr);
        }
    }
    while (curr->parent) {
        recompute_parent_vertex_mark(curr);
        curr = curr->parent;
    }
}

void HDTUFOTree::remove_edge_mark(vertex_t v) {
    auto curr = &leaves[v];
    curr->edge_mark = NONE;
    for (auto neighbor : curr->neighbors) {
        if (neighbor && neighbor->neighbors_set && neighbor->neighbors_set->contains(curr)) {
            if (!neighbor->edge_marked_siblings)
                neighbor->edge_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
            neighbor->edge_marked_siblings->insert(curr);
        }
    }
    while (curr->parent) {
        recompute_parent_edge_mark(curr);
        curr = curr->parent;
    }
}

void HDTUFOTree::recompute_parent_vertex_mark(Cluster* c) {
    assert(c->parent != nullptr);
    auto parent = c->parent;
    if (c->parent_high_fanout()) {
        auto center = (c->get_degree() == 1) ? c->get_neighbor() : c;
        parent->vertex_mark = center->get_sibling_vertex_mark();
    } else {
        parent->vertex_mark = c->get_sibling_vertex_mark();
    }
}

void HDTUFOTree::recompute_parent_edge_mark(Cluster* c) {
    assert(c->parent != nullptr);
    auto parent = c->parent;
    if (c->parent_high_fanout()) {
        auto center = (c->get_degree() == 1) ? c->get_neighbor() : c;
        parent->edge_mark = center->get_sibling_edge_mark();
    } else {
        parent->edge_mark = c->get_sibling_edge_mark();
    }
}

void HDTUFOTree::add_sibling_marks(Cluster* c1, Cluster* c2) {
    assert(c1->parent == c2->parent);
    if (c1->vertex_mark) {
        if (!c2->vertex_marked_siblings)
            c2->vertex_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
        c2->vertex_marked_siblings->insert(c1);
    }
    if (c1->edge_mark) {
        if (!c2->edge_marked_siblings)
            c2->edge_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
        c2->edge_marked_siblings->insert(c1);
    }
    if (c2->vertex_mark) {
        if (!c1->vertex_marked_siblings)
            c1->vertex_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
        c1->vertex_marked_siblings->insert(c2);
    }
    if (c2->edge_mark) {
        if (!c1->edge_marked_siblings)
            c1->edge_marked_siblings = new absl::flat_hash_set<HDTUFOCluster*>();
        c1->edge_marked_siblings->insert(c2);
    }
}

std::optional<vertex_t> HDTUFOTree::GetMarkedVertexInTree(vertex_t v) {
    vertex_t marked_vertex = leaves[v].get_root()->vertex_mark;
    if (marked_vertex == NONE) return std::nullopt;
    return marked_vertex;
}

std::optional<edge_t> HDTUFOTree::GetMarkedEdgeInTree(vertex_t v) {
    vertex_t marked_vertex = leaves[v].get_root()->edge_mark;
    if (marked_vertex == NONE) return std::nullopt;
    assert(marked_tree_edges[marked_vertex] != nullptr);
    vertex_t other_endpoint = *marked_tree_edges[marked_vertex]->begin();
    return std::pair{marked_vertex, other_endpoint};
}


vertex_t HDTUFOTree::GetSizeOfTree(vertex_t v) {
    return leaves[v].get_root()->size;
}

/* Return true if and only if there is a path from vertex u to
vertex v in the tree. */
bool HDTUFOTree::IsConnected(vertex_t u, vertex_t v) {
    return leaves[u].get_root() == leaves[v].get_root();
}
