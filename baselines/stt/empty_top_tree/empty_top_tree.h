#ifndef EMPTY_TOP_TREE_H
#define EMPTY_TOP_TREE_H
#include <stdlib.h>
#include <stdbool.h>
#include "empty_tree.h"


namespace dgbs {

typedef struct tte_node_struct tte_node;
typedef struct tte_int_node_struct tte_int_node;
typedef struct tte_leaf_node_struct tte_leaf_node;

struct tte_node_struct {
    tte_int_node *parent;
    unsigned char is_leaf:1;
    unsigned char flip:1;
    unsigned char num_boundary:2;
};

struct tte_leaf_node_struct {
    tte_node info;
    struct empty_edge *edge;
};

struct tte_int_node_struct {
    tte_node info;
    tte_node *children[2];
};

// Remove the given edge from the tree.
void tte_cut(struct empty_edge *edge);

// Add the given edge to the tree.
tte_node *tte_link(struct empty_vertex *u, struct empty_vertex *v);

// Deexpose the given vertex.
tte_node *deexpose(struct empty_vertex *vert);

// Expose the given vertex. This is the simple implementation.
tte_node *expose(struct empty_vertex *vert);

// Expose the given vertex. This is the implementation from the appendix.
tte_node *expose2(struct empty_vertex *vert);

// Find the root node of the given top tree.
tte_node *find_root(tte_node *node);

// Used for deallocating the top tree.
void destroy_top_tree_containing_edge(struct empty_edge *edge);

}

#endif
