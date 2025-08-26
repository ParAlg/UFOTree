#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/spaa_rc_tree.h"
#include "../include/spaa_rc_tree_ternarized.h"

TEST(ParallelRCTreeSuite, test_constructor){
  ParallelRCTree<int> t(3);
  /*for(int i = 0; i < t.clusters.size(); i++){
    t.clusters[i].print();
  }*/
  ASSERT_EQ(t.clusters.size(), 3);
}

/*TEST(ParallelRCTreeSuite, test_batch_insert){
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelRCTree<int> tree(n,k);

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

        parlay::sequence<Edge> batch;
        for (auto update : updates) {
            batch.push_back(update.edge);
            if (batch.size() == k) {
                tree.batch_link(batch);
                batch.clear();
                tree.verify_tree_correctness();
            }
        }
        tree.batch_link(batch);
        tree.verify_tree_correctness();

    }
}*/
TEST(ParallelRCTreeSuite, test_path_query){
    std::vector<int> test_vals = {100, 1000};
    srand(time(NULL));
    int seed = rand();
    srand(seed);
    int num_trials = 1;
    for (int n : test_vals){
        for (int trial = 0; trial < num_trials; ++trial){
            ParallelRCTree<int> tree(n, 10, std::numeric_limits<int>::max(), [](int x, int y)
                             { return std::min(x, y); });
            vertex_t u = rand() % (n - 1), v = rand() % n;
            if (v == u)
            {
                v++;
            }
            // vertex_t u = 4, v = 7;
            if (v < u)
                std::swap(u, v);
            /*std::cout << "Seed: " << seed << "\n";
          std::cout << "u = " << u << "v = " << v << "\n";*/
            int min_edge_val = std::numeric_limits<int>::max();
            parlay::sequence<std::tuple<int,int,int> > batch;
            for (int i = 0; i < n - 1; i++) {
                int new_edge = rand() % 100;
                if (i >= u && i < v)
                    min_edge_val = std::min(min_edge_val, new_edge);
                batch.push_back({i, i+1, new_edge});
                if(batch.size() == 10){
                    tree.batch_link(batch);
                    batch.clear();
                }
            }
            tree.batch_link(batch);
            // Test return of min_edge_value.
            auto returned_query = tree.path_query(u, v);
            if (returned_query != min_edge_val) {
                //print_tree(&tree);
                std::cout << seed << "\n";
                std::cout << "u = " << u << "v = " << v << "\n";
            }
            ASSERT_EQ(returned_query, min_edge_val);
        }
    }
}
TEST(ParallelRCTreeTernarizedSuite, test_path_query){
    std::vector<int> test_vals = {100, 1000};
    srand(time(NULL));
    int seed = rand();
    srand(seed);
    int num_trials = 1;
    for (int n : test_vals){
        for (int trial = 0; trial < num_trials; ++trial){
            ParallelRCTreeTernarized<int> tree(n, 10, std::numeric_limits<int>::max(), [](int x, int y)
                             { return std::min(x, y); });
            vertex_t u = rand() % (n - 1), v = rand() % n;
            if (v == u)
            {
                v++;
            }
            // vertex_t u = 4, v = 7;
            if (v < u)
                std::swap(u, v);
            /*std::cout << "Seed: " << seed << "\n";
          std::cout << "u = " << u << "v = " << v << "\n";*/
            int min_edge_val = std::numeric_limits<int>::max();
            parlay::sequence<std::tuple<int,int,int> > batch;
            for (int i = 0; i < n - 1; i++) {
                int new_edge = rand() % 100;
                if (i >= u && i < v)
                    min_edge_val = std::min(min_edge_val, new_edge);
                batch.push_back({i, i+1, new_edge});
                if(batch.size() == 10){
                    tree.batch_link(batch);
                    batch.clear();
                }
            }
            tree.batch_link(batch);
            // Test return of min_edge_value.
            auto returned_query = tree.path_query(u, v);
            if (returned_query != min_edge_val) {
                //print_tree(&tree);
                std::cout << seed << "\n";
                std::cout << "u = " << u << "v = " << v << "\n";
            }
            ASSERT_EQ(returned_query, min_edge_val);
        }
    }
}