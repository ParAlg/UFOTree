/**
 * This file contains an implementation of Kruskal's algorithm using union-find.
 */

#include "kruskal.h"

struct union_find {
    size_t *parent, *size;
};

static struct union_find *uf_make(size_t num_vertices) {
    size_t *parent, *size;
    struct union_find *result;

    parent = malloc(num_vertices * sizeof(size_t));
    if (!parent) goto fail_alloc_parent;
    size = malloc(num_vertices * sizeof(size_t));
    if (!size) goto fail_alloc_size;
    result = malloc(sizeof(struct union_find));
    if (!result) goto fail_alloc_result;

    for (size_t i = 0; i < num_vertices; ++i) {
        size[i] = 1;
        parent[i] = i;
    }

    result->parent = parent;
    result->size = size;

    return result;

fail_alloc_result:
    free(size);
fail_alloc_size:
    free(parent);
fail_alloc_parent:
    return NULL;
}

static void uf_free(struct union_find *uf) {
    free(uf->parent);
    free(uf->size);
    free(uf);
}

static size_t uf_find(struct union_find *uf, size_t x) {
    size_t parent = uf->parent[x];
    if (parent == x) return x;

    size_t root = uf_find(uf, parent);
    uf->parent[x] = root;
    return root;
}

static void uf_union(struct union_find *uf, size_t x, size_t y) {
    x = uf_find(uf, x);
    y = uf_find(uf, y);

    if (x == y) return;

    if (uf->size[x] < uf->size[y]) {
        size_t tmp = x;
        x = y;
        y = tmp;
    }

    uf->size[x] += uf->size[y];
    uf->parent[y] = x;
}

static int graph_edge_comp(const void *edge1, const void *edge2) {
    int weight1 = ((struct graph_edge *) edge1)->weight;
    int weight2 = ((struct graph_edge *) edge2)->weight;

    if (weight1 > weight2) return 1;
    if (weight1 < weight2) return -1;
    return 0;
}

bool kruskal(struct graph *graph, int *total_weight) {
    struct union_find *uf;
    uf = uf_make(graph->num_vertices);
    if (!uf) return false;

    qsort(graph->edges, graph->num_edges, sizeof(struct graph_edge), graph_edge_comp);

    *total_weight = 0;
    for (size_t i = 0; i < graph->num_edges; ++i) {
        struct graph_edge *edge = &graph->edges[i];
        size_t left = edge->left;
        size_t right = edge->right;

        if (uf_find(uf, left) != uf_find(uf, right)) {
            *total_weight += edge->weight;
            uf_union(uf, left, right);
        }
    }

    uf_free(uf);
    return true;
}
