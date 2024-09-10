#pragma once

#include "types.h"
#include "util.h"
#include <parlay/sequence.h>


using namespace std;
using namespace parlay;

enum ClusterStatus {
    NORMAL,
    ROOT,
    DEL
};

class NeighborIterator {
public:
    virtual vertex_t next() = 0;
};

class ForestStructure {
public:
    // Insert the set of vertices in V
    virtual void insert_vertices(sequence<vertex_t>& V) = 0;
    // Remove the set of vertices in V and any adjacent edges
    virtual void delete_vertices(sequence<vertex_t>& V) = 0;
    // Insert the set of edges in E
    virtual void insert_edges(sequence<Edge>& E) = 0;
    // Remove the set of edges in E
    virtual void delete_edges(sequence<Edge>& E) = 0;
    // Return the set of endpoints of the edges in E
    virtual sequence<vertex_t> get_endpoints(sequence<Edge>& E) = 0;

    // Return the edges in E that still exist in the next level up
    virtual sequence<Edge> map_edges_to_parents(sequence<Edge>& E) = 0;
    // Return the set of parents of V
    virtual sequence<vertex_t> get_parents(sequence<vertex_t>& V) = 0;
    // Set the parents of each cluster in V to the corresponding cluster in P also adding to their child count
    virtual void set_parents (sequence<vertex_t>& V, sequence<vertex_t>& P) = 0;
    // Remove the parent of each cluster in V also subtracting from the child count of the parent
    virtual void unset_parents (sequence<vertex_t>& V) = 0;
    // For each (possibly multiple) instance of a cluster in V add one to its child count
    virtual void add_children(sequence<vertex_t>& V) = 0;
    // For each (possibly multiple) instance of a cluster in V subtract one from its child count
    virtual void subtract_children(sequence<vertex_t>& V) = 0;

    // Non batch read-only helper functions
    virtual vertex_t get_degree(vertex_t v) = 0;
    virtual std::unique_ptr<NeighborIterator> get_neighbor_iterator(vertex_t v) = 0;
    virtual vertex_t get_first_neighbor(vertex_t v) = 0;
    virtual vertex_t get_other_neighbor(vertex_t v, vertex_t x) = 0;
    virtual vertex_t get_parent(vertex_t v) = 0;
    virtual void set_parent (vertex_t v, vertex_t p) = 0;
    virtual void unset_parent(vertex_t v) = 0;
    virtual vertex_t get_child_count(vertex_t v) = 0;
    virtual bool contracts(vertex_t v) = 0;

    // Batch update helper functions
    virtual vertex_t get_partner(vertex_t v) = 0;
    virtual void set_partner(vertex_t v, vertex_t p) = 0;
    virtual bool try_set_partner_atomic(vertex_t v, vertex_t p) = 0;
    virtual void unset_partner(vertex_t v) = 0;
    virtual bool is_local_max_priority(vertex_t v) = 0;
    virtual ClusterStatus get_status(vertex_t v) = 0;
    virtual void set_status(vertex_t v, ClusterStatus s) = 0;
    virtual bool try_set_status_atomic(vertex_t v, ClusterStatus s) = 0;
    void unset_status(vertex_t v);
    virtual void mark(vertex_t v) = 0;
    virtual void unmark(vertex_t v) = 0;
    virtual bool is_marked(vertex_t v) = 0;
};
