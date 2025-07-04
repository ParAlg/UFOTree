#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "parallel_ufo_tree.h"
#include "../benchmark/dynamic_trees/benchmark.h"
#include "../benchmark/parallel/parallel_benchmark.h"


using namespace dgbs;

template<typename aug_t>
bool ParallelUFOTree<aug_t>::is_valid() {
    parlay::sequence<vertex_t> clusters = tabulate(n, [&] (size_t i) { return (vertex_t) i; });
    parlay::sequence<vertex_t> next_clusters;
    for (vertex_t cluster : clusters) { // Ensure that every pair of incident vertices are in the same component
        auto iter = forests[0].get_neighbor_iterator(cluster);
        for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) { // This ensures all connectivity is correct by transitivity
            if (!connected(cluster, neighbor)) {
                std::cerr << "INCORRECT CONNECTIVITY." << std::endl;
                return false;
            }
        }
    }
    int level = 0;
    std::unordered_map<vertex_t,vertex_t> child_counts;
    while (!clusters.empty()) {
        for (vertex_t cluster : clusters) {
            if (forests[level].get_partner(cluster) != NONE) {
                std::cerr << "PARTNER FIELD IS STILL SET FOR CLUSTER " << cluster << " AT LEVEL " << level << "." << std::endl;
                return false;
            }
            if (forests[level].get_status(cluster) != NORMAL) {
                std::cerr << "STATUS FIELD IS STILL SET FOR CLUSTER " << cluster << " AT LEVEL " << level << "." << std::endl;
                return false;
            }
            if (level > 0 && child_counts[cluster] != forests[level].get_child_count(cluster)) {
                std::cerr << "INCORRECT CHILD COUNT FOR CLUSTER " << cluster << " AT LEVEL " << level << "." << std::endl;
                std::cerr << "EXPECTED " << child_counts[cluster] << " GOT " << forests[level].get_child_count(cluster) << std::endl;
                return false;
            }
            int num_combinations = 0;
            auto iter = forests[level].get_neighbor_iterator(cluster);
            for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) { // Look for invalid combinations
                if (forests[level].get_parent(cluster) != NONE && forests[level].get_parent(cluster) == forests[level].get_parent(neighbor)) {
                    num_combinations++;
                    if ((forests[level].get_degree(cluster) >= 3 || forests[level].get_degree(neighbor) >= 3)
                    && !(forests[level].get_degree(cluster) == 1 || forests[level].get_degree(neighbor) == 1)) {
                        std::cerr << "INVALID COMBINATION." << std::endl;
                        return false;
                    }
                }
            }
            if (forests[level].get_degree(cluster) == 2 && num_combinations > 1) {
                std::cerr << "INVALID COMBINATION. DEGREE 2 CLUSTER COMBINES WITH BOTH NEIGHBORS." << std::endl;
                return false;
            }
            if (forests[level].get_degree(cluster) <= 3 && !forests[level].contracts(cluster)) { // Ensure maximality of contraction
                if (forests[level].get_degree(cluster) == 1) {
                    vertex_t neighbor = forests[level].get_first_neighbor(cluster);
                    if (forests[level].get_degree(neighbor) > 2) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL. DEG 1 MUST CONTRACT WITH DEG 3+." << std::endl;
                            return false;
                        }
                    else if (!forests[level].contracts(neighbor)) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL. DEG 1 CAN CONTRACT WITH UNCONTRACTING NEIGHBOR." << std::endl;
                            return false;
                        }
                } else if (forests[level].get_degree(cluster) == 2) {
                    auto iter = forests[level].get_neighbor_iterator(cluster);
                    for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next()) {
                        if (forests[level].get_degree(neighbor) < 3 && !forests[level].contracts(neighbor)) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL. DEG 2 CAN CONTRACT WITH UNCONTRACTING NEIGHBOR." << std::endl;
                            return false;
                        }
                    }
                }
            }
        }
        child_counts.clear();
        for (vertex_t cluster : clusters) {
            vertex_t parent = forests[level].get_parent(cluster);
            if (child_counts.find(parent) == child_counts.end())
                child_counts[parent] = 1;
            else
                child_counts[parent] += 1;
        }
        clusters = forests[level++].get_parents(clusters);
    }
    return true;
}

template<typename aug_t>
int ParallelUFOTree<aug_t>::get_height(vertex_t v) {
    return forests.size();
}

template<typename aug_t>
void ParallelUFOTree<aug_t>::print_tree() {
    std::multimap<vertex_t,vertex_t> clusters;
    std::multimap<vertex_t,vertex_t> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    for (vertex_t i = 0; i < n; i++) clusters.insert({forests[0].get_parent(i), i});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << leaf << "\t " << " Parent " << (parent==NONE ? "NONE" : std::to_string(parent)) << "\t Neighbors: ";
        auto iter = forests[0].get_neighbor_iterator(leaf);
        for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next())  std::cout << neighbor << " ";
        std::cout << std::endl;
        bool in_map = false;
        for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
        if (parent != NONE && !in_map) next_clusters.insert({forests[1].get_parent(parent), parent});
    }
    clusters.swap(next_clusters);
    next_clusters.clear();
    int level = 0;
    while (!clusters.empty()) {
        level++;
        std::cout << "======================== LEVEL " << level << " ========================" << std::endl;
        for (auto entry : clusters) {
            auto cluster = entry.second;
            auto parent = entry.first;
            std::cout << "CLUSTER: " << cluster << "\t" << " Parent: " << (parent==NONE ? "NONE" : std::to_string(parent)) << "\t Neighbors: ";
            auto iter = forests[level].get_neighbor_iterator(cluster);
            for(vertex_t neighbor = iter->next(); neighbor != NONE; neighbor = iter->next())  std::cout << neighbor << " ";
            std::cout << std::endl;
            bool in_map = false;
            for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
            if (parent != NONE && !in_map) next_clusters.insert({forests[level+1].get_parent(parent), parent});
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
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

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
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

        auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_karytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t fanout = 32; // K-ary tree
        vertex_t k = 16; // batch size
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n, rand(), fanout);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
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

        auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
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

        auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
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

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
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

        auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_karytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t fanout = 32; // K-ary tree
        vertex_t k = 16; // batch size
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelUFOTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n, rand(), fanout);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
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

        auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
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

        auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of cuts.";
            }
        }
    }
}
