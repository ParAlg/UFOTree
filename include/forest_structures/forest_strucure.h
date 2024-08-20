#pragma once

#include "types.h"
#include "util.h"
#include <parlay/sequence.h>


using namespace std;
using namespace parlay;

class ForestStructure {
public:
    // Insert the set of vertices in V
    virtual void insert_vertices(sequence<vertex_t>& V) = 0;
    // Remove the set of vertices in V
    virtual void delete_vertices(sequence<vertex_t>& V) = 0;
    // Insert the set of edges in E
    virtual void insert_edges(sequence<Edge>& E) = 0;
    // Remove the set of edges in E
    virtual void delete_edges(sequence<Edge>& E) = 0;
    // Return the set of endpoints of the edges in E
    virtual sequence<vertex_t> get_endpoints(sequence<Edge>& E) = 0;

    // Return the edges in E that still exist in the next level up
    virtual sequence<Edge> filter_edges(sequence<Edge>& E) = 0;
    // Return the set of parents of V
    virtual sequence<vertex_t> get_parents(sequence<vertex_t>& V) = 0;
    // Return the set of parents of V and the number of clusters in V that are a child of each parent
    virtual sequence<pair<vertex_t,vertex_t>> count_parents(sequence<vertex_t>& V) = 0;
    // Set the parents of each cluster in V to the corresponding cluster in P
    virtual void set_parents (sequence<vertex_t>& V, sequence<vertex_t>& P) = 0;
    // Remove the parent of each cluster in V
    virtual void unset_parents (sequence<vertex_t>& V) = 0;
    // For each (possibly multiple) instance of a cluster in V add one to its child count
    virtual void add_children(sequence<vertex_t>& V) = 0;

    // Non batch helper functions
    virtual vertex_t get_degree(vertex_t v) = 0;
    virtual sequence<vertex_t> get_neighbors(vertex_t v) = 0;
    virtual void set_parent(vertex_t v, vertex_t p) = 0;
    virtual void unset_parent(vertex_t v) = 0;
    virtual vertex_t get_parent(vertex_t v) = 0;
    virtual vertex_t get_child_count(vertex_t v) = 0;
    virtual bool contracts(vertex_t v) = 0;
};
