#pragma once
#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/spaa_rc_tree_ternarized.h"

TEST(ParallelRCTreeTernarizedSuite, test_constructor){
  ParallelRCTreeTernarized<int> t(3);
  /*for(int i = 0; i < t.clusters.size(); i++){
    t.clusters[i].print();
  }*/
  ASSERT_EQ(t.clusters.size(), 6);
}

TEST(ParallelRCTreeTernarizedSuite, test_batch_insert_linked_list){
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelRCTreeTernarized<int> tree(n,k); 

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

        parlay::sequence<std::tuple<int,int,int> > batch;
        for (auto update : updates) {
            batch.push_back({update.edge.src, update.edge.dst, 0});
            if (batch.size() == k) {
                tree.batch_link(batch);
                batch.clear();
                tree.verify_tree_correctness(); 
            }
        }
        tree.batch_link(batch);
        tree.verify_tree_correctness();
        
    }
}