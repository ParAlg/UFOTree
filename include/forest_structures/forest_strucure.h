#pragma once

#include "types.h"
#include "util.h"
#include <parlay/sequence.h>


using namespace std;
using namespace parlay;

class ForestStructure {
public:
    virtual void insert_vertices(sequence<vertex_t>& V) = 0;
    virtual void delete_vertices(sequence<vertex_t>& V) = 0;
    virtual void insert_edges(sequence<Edge>& E) = 0;
    virtual void delete_edges(sequence<Edge>& E) = 0;
    virtual sequence<bool> check_edges(sequence<Edge>& E) = 0;
    virtual vertex_t get_degree(vertex_t v) = 0;

    virtual sequence<vertex_t> get_parents(sequence<vertex_t>& V) = 0;
    virtual sequence<pair<vertex_t,vertex_t>> count_parents(sequence<vertex_t>& V) = 0;
    virtual void set_parents (sequence<vertex_t>& V, sequence<vertex_t>& P) = 0;
    virtual void unset_parents (sequence<vertex_t>& V) = 0;
    virtual void add_children(sequence<vertex_t>& V) = 0;
    virtual vertex_t get_child_count(vertex_t v) = 0;
};
