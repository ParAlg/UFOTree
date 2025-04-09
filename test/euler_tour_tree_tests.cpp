#include <gtest/gtest.h>
#include <unordered_set>
#include "../baselines/parett/dynamic_trees/euler_tour_tree/treap_ett.hpp"
#include "ufo_tree.h"

using namespace dgbs;


TEST(EulerTourTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        treap::EulerTourTree<int> tree(n);

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

TEST(EulerTourTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        treap::EulerTourTree<int> tree(n);

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

TEST(EulerTourTreeSuite, subtree_query_test) {
    vertex_t n = 256;
    treap::EulerTourTree<int> tree(n);
    treap::Node<int>::aggregate_function = [] (int x, int y) { return x+y; };

    for (int i = 0; i < n-1; i++)
        tree.Link(i,i+1);
    for (int i = 0; i < n; i++)
        tree.UpdateValue(i, 1);

    for (int i = 1; i < n-1; i++) {
        ASSERT_EQ(tree.GetSubtreeAggregate(i,i+1), i+1);
        ASSERT_EQ(tree.GetSubtreeAggregate(i,i-1), n-(i+1));
    }
}
