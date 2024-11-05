#ifndef KRUSKAL_H
#define KRUSKAL_H
#include <stdlib.h>
#include <stdbool.h>

struct graph_edge {
    int weight;
    size_t left, right;
};

struct graph {
    size_t num_vertices, num_edges;
    struct graph_edge *edges;
};

// Returns true on success, false otherwise. The total weight of the MST is
// written to the given pointer.
bool kruskal(struct graph *graph, int *total_weight);

#endif
