#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "../include/parallel_ufo_tree.h"


template<typename aug_t>
bool ParallelUFOTree<aug_t>::is_valid() {
    return true;
}

template<typename aug_t>
int ParallelUFOTree<aug_t>::get_height(vertex_t v) {
    return levels.size();
}

template<typename aug_t>
void ParallelUFOTree<aug_t>::print_tree() {
}

TEST(ParallelUFOTreeSuite, batch_incremental_linkedlist_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < n-1; i++)
            edges.push_back({ids[i],ids[i+1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
            }
        }
        tree.batch_link(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_binarytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < (n-1)/2; i++) {
            edges.push_back({ids[i],ids[2*i+1]});
            edges.push_back({ids[i],ids[2*i+2]});
        }
        if (n%2 == 0) edges.push_back({ids[(n-1)/2],ids[n-1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
            }
        }
        tree.batch_link(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_star_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < n-1; i++)
            edges.push_back({ids[0],ids[i+1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
            }
        }
        tree.batch_link(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        std::vector<int> vertex_degrees(n,0);
        while (edges.size() < n-1) {
            vertex_t u = edges.size()+1;
            vertex_t v = rand() % u;
            if (vertex_degrees[v] >= 3) continue;
            edges.push_back({ids[u],ids[v]});
            vertex_degrees[u]++;
            vertex_degrees[v]++;
        }
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
            }
        }
        tree.batch_link(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_linkedlist_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < n-1; i++)
            edges.push_back({ids[i],ids[i+1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
            }
        }
        tree.batch_link(batch, len);

        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        updates.clear();
        for (auto edge : edges) updates.push_back({DELETE,edge});
        len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_cut(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
        tree.batch_cut(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_binarytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < (n-1)/2; i++) {
            edges.push_back({ids[i],ids[2*i+1]});
            edges.push_back({ids[i],ids[2*i+2]});
        }
        if (n%2 == 0) edges.push_back({ids[(n-1)/2],ids[n-1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
            }
        }
        tree.batch_link(batch, len);

        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        updates.clear();
        for (auto edge : edges) updates.push_back({DELETE,edge});
        len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_cut(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
        tree.batch_cut(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_star_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        for (vertex_t i = 0; i < n-1; i++)
            edges.push_back({ids[0],ids[i+1]});
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
            }
        }
        tree.batch_link(batch, len);

        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        updates.clear();
        for (auto edge : edges) updates.push_back({DELETE,edge});
        len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_cut(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
        tree.batch_cut(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        std::vector<Update> updates;
        parlay::sequence<Edge> edges;
        auto seed = seeds[trial];
        srand(seed);
        parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
        ids = parlay::random_shuffle(ids, parlay::random(rand()));
        std::vector<int> vertex_degrees(n,0);
        while (edges.size() < n-1) {
            vertex_t u = edges.size()+1;
            vertex_t v = rand() % u;
            if (vertex_degrees[v] >= 3) continue;
            edges.push_back({ids[u],ids[v]});
            vertex_degrees[u]++;
            vertex_degrees[v]++;
        }
        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        for (auto edge : edges) updates.push_back({INSERT,edge});

        Edge batch[k];
        vertex_t len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_link(batch, len);
                len = 0;
            }
        }
        tree.batch_link(batch, len);

        edges = parlay::random_shuffle(edges, parlay::random(rand()));
        updates.clear();
        for (auto edge : edges) updates.push_back({DELETE,edge});
        len = 0;
        for (auto update : updates) {
            batch[len++] = update.edge;
            if (len == k) {
                tree.batch_cut(batch, len);
                len = 0;
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
        tree.batch_cut(batch, len);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
    }
}

TEST(ParallelUFOTreeSuite, batch_linkedlist_performance_test) {
    vertex_t n = 1000000;
    vertex_t k = 1;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

    std::vector<Update> updates;
    parlay::sequence<Edge> edges;
    srand(time(NULL));
    parlay::sequence<vertex_t> ids = parlay::tabulate(n, [&] (vertex_t i) { return i; });
    ids = parlay::random_shuffle(ids, parlay::random(rand()));
    for (vertex_t i = 0; i < n-1; i++)
        edges.push_back({ids[i],ids[i+1]});
    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    for (auto edge : edges) updates.push_back({INSERT,edge});

    Edge* batch = (Edge*) malloc(k*sizeof(Edge));
    vertex_t len = 0;
    for (auto update : updates) {
        batch[len++] = update.edge;
        if (len == k) {
            tree.batch_link(batch, len);
            len = 0;
        }
    }
    tree.batch_link(batch, len);

    edges = parlay::random_shuffle(edges, parlay::random(rand()));
    updates.clear();
    for (auto edge : edges) updates.push_back({DELETE,edge});
    len = 0;
    for (auto update : updates) {
        batch[len++] = update.edge;
        if (len == k) {
            tree.batch_cut(batch, len);
            len = 0;
        }
    }
    tree.batch_cut(batch, len);
    free(batch);
}
