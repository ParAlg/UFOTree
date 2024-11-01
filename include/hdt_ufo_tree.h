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
    absl::flat_hash_set<HDTUFOCluster*>* vertex_marked_children;
    absl::flat_hash_set<HDTUFOCluster*>* edge_marked_children;
    // Constructor
    HDTUFOCluster();
    void insert_vertex_marked_child(HDTUFOCluster* c);
    void insert_edge_marked_child(HDTUFOCluster* c);
    void remove_vertex_marked_child(HDTUFOCluster* c);
    void remove_edge_marked_child(HDTUFOCluster* c);
    bool contracts();
    HDTUFOCluster* get_neighbor();
    HDTUFOCluster* get_root();
    int get_degree();
    bool parent_high_fanout();
    bool contains_neighbor(HDTUFOCluster* c);
    void insert_neighbor(HDTUFOCluster* c);
    void remove_neighbor(HDTUFOCluster* c);
    size_t calculate_size();
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
    void link(vertex_t u, vertex_t v) {AddEdge({u,v});};
    void cut(vertex_t u, vertex_t v) {DeleteEdge({u,v});};
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
    void remove_ancestors(Cluster* c, int start_level = 0, int intitial_size_delta = 0);
    void recluster_tree();
    void disconnect_siblings(Cluster* c, int level);
    void insert_adjacency(Cluster* u, Cluster* v);
    void remove_adjacency(Cluster* u, Cluster* v);
    void add_vertex_mark(vertex_t v);
    void add_edge_mark(vertex_t v);
    void remove_vertex_mark(vertex_t v);
    void remove_edge_mark(vertex_t v);
    void recompute_vertex_mark(Cluster* c);
    void recompute_edge_mark(Cluster* c);
};
