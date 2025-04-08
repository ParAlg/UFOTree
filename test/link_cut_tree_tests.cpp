#include <gtest/gtest.h>
#include <unordered_set>
#include "../baselines/parett/dynamic_trees/link_cut_tree/link_cut_tree.hpp"
#include "ufo_tree.h"

using namespace dgbs;
using namespace link_cut_tree;


TEST(LinkCutTreeIntSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (int j = 0; j <= i; ++j) ASSERT_TRUE(tree.connected(j,j+1));
    }
}

TEST(LinkCutTreeIntSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        for (int j = 0; j<= i; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
    }
    if (n%2 == 0) {
        tree.link((n-1)/2,n-1);
        for (int j = 0; j < (n-1)/2; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
        ASSERT_TRUE(tree.connected((n-1)/2,n-1));
    }
}

TEST(LinkCutTreeIntSuite, incremental_star_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
        for (int j = 1; j <= i+1; ++j) ASSERT_TRUE(tree.connected(0,j));
    }
}

TEST(LinkCutTreeIntSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        LinkCutTreeInt tree(n);

        auto seed = seeds[trial];
        srand(seed);
        std::vector<Edge> edges;
        while (edges.size() < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges.push_back({u,v});
                for (auto e : edges)
                    ASSERT_TRUE(tree.connected(e.src, e.dst));
            }
        }
    }
}

TEST(LinkCutTreeIntSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_FALSE(tree.connected(i,i+1));
        for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(j,j+1));
    }
}

TEST(LinkCutTreeIntSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_FALSE(tree.connected(i,2*i+1));
        ASSERT_FALSE(tree.connected(i,2*i+2));
        for (int j = i+1; j < (n-1)/2; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_FALSE(tree.connected((n-1)/2,n-1));
}

TEST(LinkCutTreeIntSuite, decremental_star_correctness_test) {
    vertex_t n = 256;
    LinkCutTreeInt tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(0,i+1);
        ASSERT_FALSE(tree.connected(0,i+1));
        for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(0,j+1));
    }
}

TEST(LinkCutTreeIntSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        LinkCutTreeInt tree(n);

        auto seed = seeds[trial];
        srand(seed);
        std::vector<Edge> edges;
        while (edges.size() < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges.push_back({u,v});
            }
        }
        for (int i = 0; i < n-1; i++) {
            auto e = edges[i];
            tree.cut(e.src, e.dst);
            ASSERT_FALSE(tree.connected(e.src,e.dst));
            for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(edges[j].src,edges[j].dst));
        }
    }
}

TEST(LinkCutTreeIntSuite, testLinkedListMaxQueries) {
    std::vector<int> test_vals = {10, 100, 1000};
    srand(time(NULL));
    int seed = rand();
    srand(seed);
    int num_trials = 100;
    for(int n : test_vals) {
        for(int trial = 0; trial < num_trials; ++trial) {
            LinkCutTreeInt tree(n);
            vertex_t u = rand() % (n-1), v = rand() % n;
            if(v == u) v++;
            if(v < u) std::swap(u,v);
            int max_edge_val = std::numeric_limits<int>::min();
            for(int i = 0; i < n-1; i++) {
                int weight = rand() % 100;
                tree.link(i, i+1, weight);
                if (i >= u && i < v) max_edge_val = std::max(max_edge_val, weight);
            }
            auto returned_query = tree.path_query(u,v).first;
            ASSERT_EQ(returned_query, max_edge_val);
        }
    }
}

TEST(LinkCutTreeIntSuite, BinaryTreeMaxQueryTest) {
    std::vector<int> test_vals = {7, 31, 127, 1023};
    srand(time(NULL));
    int seed = rand();
    seed = 0;
    srand(seed);
    int num_trials = 100;
    for(int n : test_vals) {
        for(int trial = 0; trial < num_trials; ++trial){
            LinkCutTreeInt tree(n);
            std::vector<int> path;
            for (int i = 0; i < n - 1; i = (2*i) + 2) path.push_back(i);
            vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();
            if(lower == upper) lower++;
            if(lower < upper) std::swap(lower, upper);
            auto u = path[upper], v = path[lower];
            int j = 0;
            int max_edge_val = std::numeric_limits<int>::min();
            for(int i = 0; i < (n/2); i++) { 
                auto weight1 = rand() % 100;
                auto weight2 = rand() % 100; 
                tree.link(i, (2*i) + 1, weight1);
                tree.link(i, (2*i) + 2, weight2);
                if (i == path[j]) {
                    if (i >= u && i < v) max_edge_val = std::max(max_edge_val, weight2);
                    j++;
                }
            }
            auto returned_query = tree.path_query(u,v).first;
            ASSERT_EQ(returned_query, max_edge_val);
        }
    }
}

TEST(LinkCutTreeIntSuite, PathQuerySanityCheckTest) {
    int num_trials = 1;
    int num_queries = 1000;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType q = PATH;
        auto f = [](int x, int y){return std::max(x,y);};
        int im = INT_MIN;
        UFOTree<int, int> ufo(n, q, f, f, im, im, im, im);
        LinkCutTreeInt lct(n);

        auto seed = seeds[trial];
        srand(seed);
        std::vector<Edge> edges;
        while (edges.size() < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !ufo.connected(u,v)) {
                int weight = rand() % 1000;
                ufo.link(u,v,weight);
                lct.link(u,v,weight);
                edges.push_back({u,v});
            }
        }
        for (int query = 0; query < num_queries; query++) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v) ASSERT_EQ(ufo.path_query(u,v), lct.path_query(u,v).first) << u << " " << v;
        }
        for (int i = 0; i < n/10; i++) {
            auto e = edges[i];
            ufo.cut(e.src, e.dst);
            lct.cut(e.src, e.dst);
        }
        for (int query = 0; query < num_queries; query++) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && ufo.connected(u,v))
                ASSERT_EQ(ufo.path_query(u,v), lct.path_query(u,v).first) << u << " " << v;
        }
    }
}

// =============================== REGULAR LINK CUT TREE TESTS BELOW THIS POINT ========================================

TEST(LinkCutTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (int j = 0; j <= i; ++j) ASSERT_TRUE(tree.connected(j,j+1));
    }
}

TEST(LinkCutTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        for (int j = 0; j<= i; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
    }
    if (n%2 == 0) {
        tree.link((n-1)/2,n-1);
        for (int j = 0; j < (n-1)/2; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
        ASSERT_TRUE(tree.connected((n-1)/2,n-1));
    }
}

TEST(LinkCutTreeSuite, incremental_star_correctness_test) {
    vertex_t n = 256;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
        for (int j = 1; j <= i+1; ++j) ASSERT_TRUE(tree.connected(0,j));
    }
}

TEST(LinkCutTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        LinkCutTree tree(n);

        auto seed = seeds[trial];
        srand(seed);
        std::vector<Edge> edges;
        while (edges.size() < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges.push_back({u,v});
                for (auto e : edges) ASSERT_TRUE(tree.connected(e.src, e.dst));
            }
        }
    }
}

TEST(LinkCutTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 128;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_FALSE(tree.connected(i,i+1));
        for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(j,j+1));
    }
}

TEST(LinkCutTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_FALSE(tree.connected(i,2*i+1));
        ASSERT_FALSE(tree.connected(i,2*i+2));
        for (int j = i+1; j < (n-1)/2; ++j) ASSERT_TRUE(tree.connected(j,2*j+1) && tree.connected(j,2*j+2));
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_FALSE(tree.connected((n-1)/2,n-1));
}

TEST(LinkCutTreeSuite, decremental_star_correctness_test) {
    vertex_t n = 256;
    LinkCutTree tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(0,i+1);
        ASSERT_FALSE(tree.connected(0,i+1));
        for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(0,j+1));
    }
}

TEST(LinkCutTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        LinkCutTree tree(n);

        auto seed = seeds[trial];
        srand(seed);
        std::vector<Edge> edges;
        while (edges.size() < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges.push_back({u,v});
            }
        }
        for (int i = 0; i < n-1; i++) {
            auto e = edges[i];
            tree.cut(e.src, e.dst);
            ASSERT_FALSE(tree.connected(e.src,e.dst));
            for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(edges[j].src,edges[j].dst));
        }
    }
}
