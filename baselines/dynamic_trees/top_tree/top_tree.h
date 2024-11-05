#ifndef TOP_TREE_H
#define TOP_TREE_H
#include <stdlib.h>
#include <stdbool.h>
#include "tree.h"

typedef struct tt_node_struct tt_node;
typedef struct tt_int_node_struct tt_int_node;
typedef struct tt_leaf_node_struct tt_leaf_node;

struct tt_node_struct {
    tt_int_node *parent;
    int spine_weight;
    unsigned char is_leaf:1;
    unsigned char flip:1;
    unsigned char num_boundary:2;
};

struct tt_leaf_node_struct {
    tt_node info;
    struct edge *edge;
};

struct tt_int_node_struct {
    tt_node info;
    tt_node *children[2];
};

// Remove the given edge from the tree.
void tt_cut(struct edge *edge);

// Add the given edge to the tree.
tt_node *tt_link(struct vertex *u, struct vertex *v, int weight);

// Deexpose the given vertex.
tt_node *deexpose(struct vertex *vert);

// Expose the given vertex. This is the simple implementation.
tt_node *expose(struct vertex *vert);

// Expose the given vertex. This is the implementation from the appendix.
tt_node *expose2(struct vertex *vert);

// Finds an edge with maximum weight between the two exposed vertices, given the
// root of a top tree.
tt_leaf_node *find_maximum(tt_node *root);

// Find the root node of the given top tree.
tt_node *find_root(tt_node *node);

// Used for deallocating the top tree.
void destroy_top_tree_containing_edge(struct edge *edge);

#endif
