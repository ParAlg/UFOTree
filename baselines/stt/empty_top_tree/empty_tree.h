#ifndef EMPTY_TREE_H
#define EMPTY_TREE_H
#include <stdlib.h>
#include <stdbool.h>
#include <absl/container/flat_hash_map.h>
#include "../../../util/util.h"


namespace dgbs {

static absl::flat_hash_map<std::pair<struct empty_vertex*, struct empty_vertex*>, struct empty_edge*> empty_edges;
static long long tte_changes = 0;

struct empty_tree {
    size_t num_vertices;
    struct empty_vertex* vertices;
};

struct empty_vertex {
    struct empty_edge *first_edge;
    bool is_exposed;
};

// Each edge is part of two linked lists. One for `left` and one for `right`.
struct empty_edge {
    void *user_data;
    struct empty_vertex *endpoints[2];
    struct empty_edge *prev[2];
    struct empty_edge *next[2];
};

struct empty_tree create_empty_tree(size_t num_vertices);
void add_edge(struct empty_edge *allocation, struct empty_vertex *left, struct empty_vertex *right);
void destroy_edge(struct empty_edge *edge);
void destroy_tree(struct empty_tree *to_destroy);
bool has_at_most_one_incident_edge(struct empty_vertex *vertex);

}

#endif
