#include "top_tree.h"
#include "tree.h"
#include <stdio.h>

static int count_leaves_with(tt_node *node, struct vertex *vert) {
    if (node->is_leaf) {
        tt_leaf_node *leaf = (tt_leaf_node *) node;
        return leaf->edge->endpoints[0] == vert || leaf->edge->endpoints[1] == vert;
    } else {
        tt_int_node *int_node = (tt_int_node *) node;
        return count_leaves_with(int_node->children[0], vert)
            + count_leaves_with(int_node->children[1], vert);
    }
}

static size_t get_degree(struct vertex *vert) {
    size_t count = 0;
    struct edge *edge = vert->first_edge;
    while (edge) {
        count += 1;
        int j = edge->endpoints[1] == vert;
        edge = edge->next[j];
    }
    return count;
}

struct cii_res {
    struct vertex *left, *mid, *right;
};
struct cii_res check_invariants_inner(tt_node *node) {
    struct cii_res res;
    res.left = NULL;
    res.mid = NULL;
    res.right = NULL;
    if (node->is_leaf) {
        tt_leaf_node *leaf = (tt_leaf_node *) node;
        struct vertex *ep_left = leaf->edge->endpoints[leaf->info.flip];
        struct vertex *ep_right = leaf->edge->endpoints[!leaf->info.flip];
        if (ep_left->is_exposed || !has_at_most_one_incident_edge(ep_left)) {
            res.left = ep_left;
        }
        if (ep_right->is_exposed || !has_at_most_one_incident_edge(ep_right)) {
            res.right = ep_right;
        }

        if ((!!res.left) + (!!res.mid) + (!!res.right) != node->num_boundary) {
            printf("num_boundary mismatch %d\n", node->is_leaf);
        }
    } else {
        tt_int_node *int_node = (tt_int_node *) node;
        struct cii_res bl = check_invariants_inner(int_node->children[0]);
        struct cii_res br = check_invariants_inner(int_node->children[1]);

        struct vertex *bl_rightmost = bl.right ? bl.right : bl.mid;
        struct vertex *br_rightmost = br.right ? br.right : br.mid;
        struct vertex *bl_leftmost = bl.left ? bl.left : bl.mid;
        struct vertex *br_leftmost = br.left ? br.left : br.mid;
        if (bl_rightmost != br_leftmost) {
            printf("left->rightmost != right->leftmost\n");
        }

        int leaves = count_leaves_with(node, bl_rightmost);
        int deg = get_degree(bl_rightmost);
        if (bl_rightmost->is_exposed || (leaves < deg)) {
            res.mid = bl_rightmost;
        }
        if (bl_leftmost != bl_rightmost) {
            res.left = bl_leftmost;
        }
        if (br_leftmost != br_rightmost) {
            res.right = br_rightmost;
        }

        if (int_node->info.flip) {
            struct vertex *tmp = res.left;
            res.left = res.right;
            res.right = tmp;
        }

        if ((!!res.left) + (!!res.mid) + (!!res.right) != node->num_boundary) {
            printf("num_boundary mismatch %d\n", node->is_leaf);
        }
    }

    if (!node->parent) {
        if ((res.left && !res.left->is_exposed)
                || (res.mid && !res.mid->is_exposed)
                || (res.right && !res.right->is_exposed)) {
            printf("root is wrong\n");
        }
    }

    return res;
}

void check_invariants(tt_node *node) {
    if (node) check_invariants_inner(find_root(node));
}
