#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "top_tree.h"
#include "tree.h"
#include "kruskal.h"
#include "debug.h"

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

static void get_neighbours(struct vertex *vertices, size_t i, size_t *buf) {
    struct edge *edge = vertices[i].first_edge;
    while (edge) {
        int j = edge->endpoints[1] == vertices+i;
        *buf++ = edge->endpoints[!j] - vertices;
        edge = edge->next[j];
    }
}

static int sum_over_adjacent_edges(struct vertex *vertices, size_t i) {
    int sum = 0;
    struct edge *edge = vertices[i].first_edge;
    while (edge) {
        sum += edge->weight;
        int j = edge->endpoints[1] == vertices+i;
        edge = edge->next[j];
    }
    return sum;
}

static int compare_size_t(const void* a, const void* b) {
     size_t int_a = * ( (size_t*) a );
     size_t int_b = * ( (size_t*) b );

     if ( int_a == int_b ) return 0;
     else if ( int_a < int_b ) return -1;
     else return 1;
}

static int run_command_mode() {
    char *cmd_buf = malloc(sizeof(char) * 101);
    size_t *print_buf = NULL;
    size_t print_buf_len = 0;
    if (!cmd_buf) return 1;

    int N;
    int num_match = scanf(" %d", &N);
    if (num_match < 1) {
        free(cmd_buf);
        return 1;
    }

    struct tree T;
    T = create_tree(N);
    if (!T.vertices) {
        free(cmd_buf);
        return 1;
    }

    struct vertex *v = T.vertices;
    struct vertex *end = v + T.num_vertices;

    int ret = 0;
    while (true) {
        num_match = scanf(" %100s", cmd_buf);
        if (num_match < 1) break;

        if (strcmp(cmd_buf, "link") == 0) {
            int ui, vi, new_weight;
            num_match = scanf(" %d %d %d", &ui, &vi, &new_weight);
            if (num_match < 3) {
                ret = 1;
                goto exit;
            }
            if (ui >= T.num_vertices || vi >= T.num_vertices) {
                printf("Index too large.\n");
                ret = 1;
                goto exit;
            }

            tt_node *root1 = expose(&v[ui]);
            tt_node *root2 = expose(&v[vi]);

            tt_leaf_node *max_edge = NULL;
            bool insert_link = false;
            if (root1 && (root1 == root2)) {
                if (new_weight < root1->spine_weight) {
                    insert_link = true;
                    max_edge = find_maximum(root1);
                } else {
                    printf("Skipping link %d-%d\n", ui, vi);
                }
            } else {
                insert_link = true;
            }
            deexpose(&v[ui]);
            deexpose(&v[vi]);

            if (max_edge) {
                struct vertex *ep0 = max_edge->edge->endpoints[0];
                struct vertex *ep1 = max_edge->edge->endpoints[1];
                printf("Cutting link %d-%d\n", (int) (ep0 - v), (int) (ep1 - v));
                tt_cut(max_edge->edge);
            }

            if (insert_link) {
                printf("Inserting link %d-%d\n", ui, vi);
                tt_node *new_root = tt_link(&v[ui], &v[vi], new_weight);
                if (!new_root) {
                    printf("Allocation failure.\n");
                    ret = 1;
                    goto exit;
                }
            }
        } else if (strcmp(cmd_buf, "query") == 0) {
            int ui, vi;
            num_match = scanf(" %d %d", &ui, &vi);
            if (num_match < 2) {
                ret = 1;
                goto exit;
            }
            if (ui >= T.num_vertices || vi >= T.num_vertices) {
                printf("Index too large.\n");
                ret = 1;
                goto exit;
            }

            tt_node *root1 = expose(&v[ui]);
            tt_node *root2 = expose(&v[vi]);

            if (root1 && root1 == root2) {
                tt_leaf_node *best = find_maximum(root1);
                struct edge *edge = best->edge;

                bool match_u = edge->endpoints[0] == &v[ui] || edge->endpoints[1] == &v[ui];
                bool match_v = edge->endpoints[0] == &v[vi] || edge->endpoints[1] == &v[vi];
                if (match_u && match_v) {
                    printf("%d-%d = in MST\n", ui, vi);
                } else {
                    printf("%d-%d = not in MST\n", ui, vi);
                }
            } else {
                printf("%d-%d = not connected\n", ui, vi);
            }

            deexpose(&v[ui]);
            deexpose(&v[vi]);
        } else if (strcmp(cmd_buf, "print") == 0) {
            for (int i = 0; i < T.num_vertices; ++i) {
                if (v[i].first_edge) {
                    printf("%d: ", i);
                    size_t deg = get_degree(&v[i]);
                    if (print_buf_len < deg) {
                        size_t new_len = deg;
                        if (new_len < 2*print_buf_len) new_len = 2*print_buf_len;
                        free(print_buf);
                        print_buf_len = new_len;
                        print_buf = malloc(sizeof(size_t) * new_len);
                        if (!print_buf) {
                            print_buf_len = 0;
                            ret = 1;
                            goto exit;
                        }
                    }
                    get_neighbours(v, i, print_buf);
                    qsort(print_buf, deg, sizeof(size_t), compare_size_t);
                    printf("%d", (int)print_buf[0]);
                    for (int j = 1; j < deg; ++j)
                        printf(", %d", (int)print_buf[j]);
                    printf("\n");
                }
            }
        } else {
            printf("Invalid command %s.\n", cmd_buf);
            ret = 1;
            goto exit;
        }
    }

exit:
    for (struct vertex *vert = v; vert < end; ++vert)
        destroy_top_tree_containing_edge(vert->first_edge);
    destroy_tree(&T);
    free(cmd_buf);
    free(print_buf);

    return ret;
}

static int run_compare_mode(size_t num_vertices, size_t num_edges) {
    int ret = 1;
    size_t *missing_edges;
    struct graph g;
    struct tree t;

    missing_edges = malloc(num_vertices * (num_vertices - 1) * sizeof(size_t));
    if (!missing_edges) goto fail_alloc_missing_edges;

    g.num_vertices = num_vertices;
    g.num_edges = num_edges;
    g.edges = malloc(num_edges * sizeof(struct graph_edge));
    if (!g.edges) goto fail_alloc_graph;

    t = create_tree(num_vertices);
    if (!t.vertices) goto fail_alloc_tree;

    printf("Generating random graph with %d vertices and %d edges\n", (int)num_vertices, (int)num_edges);

    size_t max_edges = (num_vertices * (num_vertices - 1)) / 2;
    if (num_edges > max_edges) {
        printf("Too many edges.\n");
        goto fail_max_edges;
    }

    // Generate all possible edges.
    size_t *next_missing_edge = missing_edges;
    for (size_t i = 0; i < num_vertices; ++i) {
        for (size_t j = i+1; j < num_vertices; ++j) {
            *(next_missing_edge++) = i;
            *(next_missing_edge++) = j;
        }
    }

    // Pick `num_edges` from full list without replacement.
    for (size_t i = 0; i < num_edges; ++i) {
        size_t remaining_edges = max_edges - i;
        int next_edge = rand() % remaining_edges;

        g.edges[i].left = missing_edges[2*next_edge];
        g.edges[i].right = missing_edges[2*next_edge+1];
        g.edges[i].weight = rand() % 100;

        missing_edges[2*next_edge] = missing_edges[2*remaining_edges-2];
        missing_edges[2*next_edge+1] = missing_edges[2*remaining_edges-1];
    }

    free(missing_edges);
    missing_edges = NULL;

    // Generate MST using a top tree.
    printf("Graph generated. Running algorithms.\n");
    for (size_t i = 0; i < num_edges; ++i) {
        size_t a = g.edges[i].left;
        size_t b = g.edges[i].right;
        int new_weight = g.edges[i].weight;

        tt_node *root1 = expose(&t.vertices[a]);
        tt_node *root2 = expose(&t.vertices[b]);

        tt_leaf_node *max_edge = NULL;
        bool insert_link = false;
        if (root1 && (root1 == root2)) {
            if (new_weight < root1->spine_weight) {
                insert_link = true;
                max_edge = find_maximum(root1);
            }
        } else {
            insert_link = true;
        }
        deexpose(&t.vertices[a]);
        deexpose(&t.vertices[b]);
        check_invariants(root1);
        check_invariants(root2);

        if (max_edge) tt_cut(max_edge->edge);
        if (insert_link) {
            tt_node *new_root = tt_link(&t.vertices[a], &t.vertices[b], new_weight);
            if (!new_root) {
                printf("Allocation failure.\n");
                goto fail_top_tree;
            }
            check_invariants(new_root);
        }
    }

    // Compute weight by inspecting tree.
    int top_tree_total_weight = 0;
    for (int i = 0; i < t.num_vertices; ++i) {
        top_tree_total_weight += sum_over_adjacent_edges(t.vertices, i);
    }
    top_tree_total_weight = top_tree_total_weight / 2;

    // Run Kruskal's algorithm.
    int kruskal_total_weight = 0;
    if (!kruskal(&g, &kruskal_total_weight)) {
        printf("Allocation failure.\n");
        goto fail_kruskal;
    }

    // Print result.
    if (top_tree_total_weight == kruskal_total_weight) {
        printf("Both algorithms agree!\nThe MST has weight %d.\n", top_tree_total_weight);
        ret = 0;
    } else {
        printf("The algorithms disagree.\nTop tree says %d.\nKruskal says %d.\n",
                top_tree_total_weight,
                kruskal_total_weight);
        ret = 1;
    }

fail_kruskal:
fail_top_tree:
    for (struct vertex *vert = t.vertices; vert < t.vertices + num_vertices; ++vert)
        destroy_top_tree_containing_edge(vert->first_edge);
fail_max_edges:
    destroy_tree(&t);
fail_alloc_tree:
    free(g.edges);
fail_alloc_graph:
    free(missing_edges);
fail_alloc_missing_edges:
    return ret;
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        printf("Commandline arguments missing.\n");
        return 1;
    }

    if (strcmp(argv[1], "command") == 0) {
        return run_command_mode();
    }

    if (strcmp(argv[1], "compare") == 0) {
        if (argc <= 3) {
            printf("Too few arguments. Must specify number of vertices and edges in random graph.\n");
            return 1;
        }

        int num_vertices = atoi(argv[2]);
        int num_edges = atoi(argv[3]);
        int repeats = 1;

        if (argc >= 5) {
            repeats = atoi(argv[4]);
        }

        srand(time(NULL));
        for (int i = 0; i < repeats; ++i) {
            int res = run_compare_mode(num_vertices, num_edges);
            if (res != 0) return res;
        }
        return 0;
    }

    printf("Unknown mode `%s`. Must be `command` or `compare`.", argv[1]);
    return 1;
}
