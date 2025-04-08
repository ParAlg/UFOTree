#include "empty_top_tree.h"
#include "empty_tree.h"
#include <limits.h>

#define is_point_e(n) (((tte_node *) n)->num_boundary < 2)
#define is_path_e(n) (((tte_node *) n)->num_boundary == 2)


namespace dgbs {

// static inline int int_max(int a, int b) {
//     if (a < b)
//         return b;
//     else
//         return a;
// }
// Replacement for Fail label in original code that was erroring.
static void fail(tte_leaf_node *T_edge, tte_int_node *Tu_new, tte_int_node *Tv_new){
    free(T_edge); free(Tu_new); free(Tv_new);
}

tte_node *find_root(tte_node *node) {
    if (node)
        while (node->parent)
            node = (tte_node*) node->parent;
    return node;
}

static void push_flip(tte_int_node *node) {
    if (node->info.flip) {
        node->info.flip = false;

        tte_node* tmp = node->children[0];
        node->children[0] = node->children[1];
        node->children[1] = tmp;
        tte_changes += 2;
        node->children[0]->flip ^= true;
        node->children[1]->flip ^= true;
    }
}

static bool has_left_boundary(tte_node *node) {
    if (node->is_leaf) {
        tte_leaf_node *leaf_node = (tte_leaf_node *) node;
        struct empty_vertex *endp = leaf_node->edge->endpoints[node->flip];
        return endp->is_exposed || !has_at_most_one_incident_edge(endp);
    } else {
        tte_int_node *int_node = (tte_int_node *) node;
        tte_node *child = int_node->children[node->flip];
        return is_path_e(child);
    }
}

static bool has_right_boundary(tte_node *node) {
    if (node->is_leaf) {
        tte_leaf_node *leaf_node = (tte_leaf_node *) node;
        struct empty_vertex *endp = leaf_node->edge->endpoints[!node->flip];
        return endp->is_exposed || !has_at_most_one_incident_edge(endp);
    } else {
        tte_int_node *int_node = (tte_int_node *) node;
        tte_node *child = int_node->children[!node->flip];
        return is_path_e(child);
    }
}

static bool has_middle_boundary(tte_node *node) {
    if (node->is_leaf || node->num_boundary == 0)
        return false;
    tte_int_node *int_node = (tte_int_node *) node;
    bool left_path = is_path_e(int_node->children[0]);
    bool right_path = is_path_e(int_node->children[1]);

    return node->num_boundary - left_path - right_path;
}

static tte_node *get_sibling(tte_node *node) {
    tte_int_node *parent = node->parent;
    if (!parent) return NULL;
    int j = parent->children[0] == node;
    return parent->children[j];
}

static void rotate_up(tte_node *node) {
    tte_int_node *parent = node->parent;
    tte_int_node *grandparent = parent->info.parent;
    tte_node *grandparent_node = (tte_node *) grandparent;
    tte_node *sibling = get_sibling(node);
    tte_node *uncle = get_sibling((tte_node *) parent);

    push_flip(grandparent);
    push_flip(parent);

    bool uncle_is_left_child = grandparent->children[0] == uncle;
    bool sibling_is_left_child = parent->children[0] == sibling;
    bool to_same_sides = uncle_is_left_child == sibling_is_left_child;
    bool sibling_is_path_e = is_path_e(sibling);
    bool uncle_is_path_e = is_path_e(uncle);
    bool gp_is_path_e = is_path_e(grandparent);
    bool new_parent_is_path_e, flip_new_parent, flip_grandparent;

    if (to_same_sides && sibling_is_path_e) {
        // Rotation on a path.
        bool gp_middle = has_middle_boundary(grandparent_node);
        new_parent_is_path_e = gp_middle || uncle_is_path_e;
        flip_new_parent = false;
        flip_grandparent = false;
        if (gp_middle && !gp_is_path_e) {
            tte_int_node *ggp = grandparent_node->parent;
            if (ggp) {
                bool gp_is_left_child = ggp->children[0] == grandparent_node;
                flip_grandparent = gp_is_left_child == uncle_is_left_child;
            }
        }
    } else {
        // Rotation on a star.
        if (!to_same_sides) {
            new_parent_is_path_e = sibling_is_path_e || uncle_is_path_e;
            flip_new_parent = sibling_is_path_e;
            flip_grandparent = sibling_is_path_e;
            node->flip ^= true;
        } else {
            new_parent_is_path_e = uncle_is_path_e;
            flip_new_parent = false;
            flip_grandparent = false;
            sibling->flip ^= true;
        }
    }

    parent->children[uncle_is_left_child] = sibling;
    parent->children[!uncle_is_left_child] = uncle;
    parent->info.flip = flip_new_parent;
    parent->info.num_boundary = new_parent_is_path_e + 1;
    tte_changes += 2;

    grandparent->children[uncle_is_left_child] = node;
    grandparent->children[!uncle_is_left_child] = (tte_node *) parent;
    grandparent->info.flip = flip_grandparent;
    tte_changes += 2;

    node->parent = grandparent;
    uncle->parent = parent;
    tte_changes += 2;
}

static tte_node *splay_step(tte_node *node) {
    // Uses loop rather than recursion.
    while (true) {
        tte_int_node *parent = node->parent;
        if (!parent) return NULL;
        tte_int_node *gparent = parent->info.parent;
        if (!gparent) return NULL;

        if (is_point_e(node) && is_point_e(gparent)) {
            rotate_up(node);
            return (tte_node *) gparent;
        }

        tte_int_node *ggparent = gparent->info.parent;
        if (!ggparent) return NULL;
        if (is_path_e(parent) && (is_path_e(gparent) || is_point_e(ggparent))) {
            push_flip(gparent);
            push_flip(parent);
            bool node_is_left = parent->children[0] == node;
            bool parent_is_left = gparent->children[0] == (tte_node*) parent;
            bool gparent_is_left = ggparent->children[0] == (tte_node*) gparent;
            if (node_is_left == parent_is_left) {
                rotate_up(node);
                return (tte_node*) gparent;
            }
            if (parent_is_left == gparent_is_left) {
                rotate_up((tte_node*) parent);
                return (tte_node*) ggparent;
            }
            // At this point, `node_is_left == gparent_is_left` is true
            rotate_up(get_sibling(node));
            rotate_up((tte_node*) parent);
            return (tte_node*) ggparent;
        }

        node = (tte_node*) parent;
    }
}

static void semi_splay(tte_node *node) {
    while (!node)
        node = splay_step(node);
}
static void full_splay(tte_node *node) {
    while (true) {
        tte_node *top = splay_step(node);
        if (!top) return;
        splay_step(top);
    }
}

static tte_node *find_consuming_node(struct empty_vertex *vert) {
    struct empty_edge *start = vert->first_edge;

    if (!start) return NULL;
    tte_node *node = (tte_node*) start->user_data;
    semi_splay(node);
    if (has_at_most_one_incident_edge(vert)) return node;

    bool is_left = (start->endpoints[0] == vert) != node->flip;
    bool is_middle = false;
    bool is_right = (start->endpoints[1] == vert) != node->flip;
    tte_node *last_middle_node = NULL;

    while (node->parent) {
        tte_int_node *parent = node->parent;
        bool is_left_child = parent->children[0] == node;

        // Compute where v is in the parent, taking the parent's
        // flip into account.
        is_middle = is_left_child
            ? (is_right || (is_middle && !has_right_boundary(node)))
            : (is_left || (is_middle && !has_left_boundary(node)));
        is_left = (is_left_child != parent->info.flip) && !is_middle;
        is_right = (is_left_child == parent->info.flip) && !is_middle;

        // Go up to the parent
        node = (tte_node*) parent;

        // If v is in the middle, then it could be the consuming node.
        if (is_middle) {
            if (!has_middle_boundary(node)) {
                // This only happens if the vertex is not exposed.
                return node;
            }
            last_middle_node = node;
        }
    }
    // This only happens when the vertex is exposed.
    return last_middle_node;
}

static tte_node *prepare_expose(tte_node *consuming_node) {
    tte_node *node = consuming_node;
    while (node->parent) {
        tte_int_node *parent = node->parent;
        if (is_point_e(node)) {
            node = (tte_node*) parent;
        } else {
            tte_int_node *int_node = (tte_int_node*) node;

            push_flip(parent);
            push_flip(int_node);

            tte_node *sibling = get_sibling(node);
            int sibling_idx = parent->children[1] == sibling;
            tte_node *same_side_child = int_node->children[sibling_idx];
            if (is_path_e(same_side_child) || is_point_e(sibling)) {
                // Case (a),(b),(c),(d)
                tte_node *other_side_child = int_node->children[1-sibling_idx];
                rotate_up(other_side_child);
                if (node == consuming_node) {
                    // Case (a),(b)
                    consuming_node = (tte_node*) parent;
                }
                node = (tte_node*) parent;
            } else {
                tte_node *uncle = get_sibling((tte_node*) parent);
                tte_int_node *gparent = parent->info.parent;
                int uncle_idx = gparent->children[1] == uncle;

                if (sibling_idx == uncle_idx) {
                    // Case (e)
                    rotate_up(node);
                } else {
                    // Case (f)
                    rotate_up(sibling);
                }
            }
        }
    }
    return consuming_node;
}

tte_node *expose_prepared(tte_node *consuming_node) {
    bool from_left = false;
    bool from_right = false;
    tte_node *node = consuming_node;
    while (true) {
        node->num_boundary += 1;
        tte_int_node *parent = node->parent;
        if (!parent) return node;

        bool is_left_child_of_parent = parent->children[0] == node;
        bool is_right_child_of_parent = !is_left_child_of_parent;

        if ((is_left_child_of_parent && from_right)
                || (is_right_child_of_parent && from_left)) {
            node->flip ^= true;
        }

        from_left = is_left_child_of_parent != parent->info.flip;
        from_right = is_right_child_of_parent != parent->info.flip;
        node = (tte_node*) parent;
    }
}

tte_node *expose2(struct empty_vertex *vert) {
    tte_node *consuming_node = find_consuming_node(vert);
    if (consuming_node) {
        consuming_node = prepare_expose(consuming_node);
        tte_node *root = expose_prepared(consuming_node);
        vert->is_exposed = true;
        return root;
    } else {
        // The vertex has degree zero.
        vert->is_exposed = true;
        return NULL;
    }
}

tte_node *expose(struct empty_vertex *vert) {
    tte_node *node = find_consuming_node(vert); // this contains a semi_splay
    if (!node) {
        // The vertex has degree zero.
        vert->is_exposed = true;
        return NULL;
    }

    while (is_path_e(node)) { // rotate_up until consuming node is a point cluster
        tte_int_node *int_node = (tte_int_node*) node;
        tte_int_node *parent = node->parent;
        push_flip(int_node);
        int node_idx = parent->children[1] == node;
        rotate_up(int_node->children[node_idx]);
        node = (tte_node*) parent;
    }

    full_splay(node);

    // Now depth(node)<=1, and node is the consuming point cluster.
    tte_node *root = NULL;
    while (node) {
        root = node;
        root->num_boundary += 1;
        node = (tte_node*) root->parent;
    }
    vert->is_exposed = true;
    return root;
}

tte_node *deexpose(struct empty_vertex *vert) {
    tte_node *root = NULL;
    tte_node *node = find_consuming_node(vert);
    while (node) {
        root = node;
        root->num_boundary -= 1;
        node = (tte_node*) root->parent;
    }
    vert->is_exposed = false;
    return root;
}

tte_node *tte_link(struct empty_vertex *u, struct empty_vertex *v) {
    tte_leaf_node *T_edge = NULL;
    tte_int_node *Tu_new = NULL;
    tte_int_node *Tv_new = NULL;
    struct empty_edge *edge = NULL;

    T_edge = (tte_leaf_node*) malloc(sizeof(tte_leaf_node));
    if (!T_edge) return NULL;
    if (u->first_edge) {
        Tu_new = (tte_int_node*) malloc(sizeof(tte_int_node));
        if (!Tu_new) fail(T_edge, Tu_new, Tv_new);
    }
    if (v->first_edge) {
        Tv_new = (tte_int_node*) malloc(sizeof(tte_int_node));
        if (!Tv_new) fail(T_edge, Tu_new, Tv_new);
    }
    edge = (struct empty_edge*) malloc(sizeof(struct empty_edge));
    if (!edge) fail(T_edge, Tu_new, Tv_new);

    tte_node *Tu = expose(u);
    if (Tu && has_left_boundary(Tu))
        Tu->flip = !Tu->flip;
    u->is_exposed = false;

    tte_node *Tv = expose(v);
    if (Tv && has_right_boundary(Tv))
        Tv->flip = !Tv->flip;
    v->is_exposed = false;

    // Init T_edge
    add_edge(edge, u, v);
    edge->user_data = T_edge;
    T_edge->edge = edge; tte_changes++;
    T_edge->info.parent = NULL;
    T_edge->info.is_leaf = true;
    T_edge->info.flip = 0;
    T_edge->info.num_boundary = (!!Tu) + (!!Tv);
    tte_node *T = (tte_node*) T_edge;

    if (Tu) {
        Tu_new->info.parent = NULL;
        Tu_new->info.is_leaf = false;
        Tu_new->info.flip = 0;
        Tu_new->info.num_boundary = !!Tv;
        Tu_new->children[0] = Tu;
        Tu_new->children[1] = T;
        Tu->parent = Tu_new;
        T->parent = Tu_new;
        tte_changes += 4;
        T = (tte_node*) Tu_new;
    }
    if (Tv) {
        Tv_new->info.parent = NULL;
        Tv_new->info.is_leaf = false;
        Tv_new->info.flip = 0;
        Tv_new->info.num_boundary = 0;
        Tv_new->children[0] = T;
        Tv_new->children[1] = Tv;
        T->parent = Tv_new;
        Tv->parent = Tv_new;
        tte_changes += 4;
        T = (tte_node*) Tv_new;
    }

    return T;
    
}


static void delete_all_ancestors(tte_node *node) {
    tte_node *p = (tte_node*) node->parent;
    if (p) {
        tte_node *s = get_sibling(node);
        delete_all_ancestors(p);
        s->parent = NULL;
        tte_changes++;
    }
    free(node);
}

void tte_cut(struct empty_edge *edge) {
    tte_node *node = (tte_node*) edge->user_data;

    struct empty_vertex *u = edge->endpoints[0];
    struct empty_vertex *v = edge->endpoints[1];
    full_splay(node);
    // now depth(e)<=2, and if e is a leaf edge, depth(e)<=1
    delete_all_ancestors(node);
    destroy_edge(edge);
    u->is_exposed = true;
    v->is_exposed = true;

    deexpose(u);
    deexpose(v);
}

static void free_top_tree(tte_node *node) {
    if (node->is_leaf) {
        tte_leaf_node *leaf = (tte_leaf_node *) node;
        leaf->edge->user_data = NULL;
    } else {
        tte_int_node *int_node = (tte_int_node *) node;
        free_top_tree(int_node->children[0]);
        free_top_tree(int_node->children[1]);
    }
    free(node);
}
void destroy_top_tree_containing_edge(struct empty_edge *edge) {
    if (!edge) return;
    if (!edge->user_data) return;
    free_top_tree(find_root((tte_node*) edge->user_data));
}

}
