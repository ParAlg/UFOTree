#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "../util/types.h"
#include "../include/parallel_topology_tree.h"
#include "../benchmark/dynamic_trees/benchmark.h"
#include "../benchmark/parallel/parallel_benchmark.h"


using namespace dgbs;

template<typename aug_t>
bool ParallelTopologyTree<aug_t>::is_valid() {
    std::unordered_set<ParallelTopologyCluster<aug_t>*> clusters;
    std::unordered_set<ParallelTopologyCluster<aug_t>*> next_clusters;
    for (int i = 0; i < n; i++) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaves[i].neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor && leaves[i].get_root() != neighbor->get_root()) {
                std::cerr << "INCORRECT CONNECTIVITY." << std::endl;
                return false;
            }
    for (int i = 0; i < n; i++) clusters.insert(&leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor && !neighbor->contains_neighbor(cluster)) {
                    std::cerr << "INCORRECT ADJACENCY." << std::endl;
                    return false;
                }
            if (!cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts()) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                            return false;
                        }
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 3) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                            return false;
                        }
                } else if (cluster->get_degree() == 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 2) {
                            std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                            return false;
                        }
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

template<typename aug_t>
int ParallelTopologyTree<aug_t>::get_height(vertex_t v) {
    int height = 0;
    ParallelTopologyCluster<aug_t>* curr = &leaves[v];
    while (curr) {
        height++;
        curr = curr->parent;
    }
    return height;
}

template<typename aug_t>
void ParallelTopologyTree<aug_t>::print_tree() {
    std::multimap<ParallelTopologyCluster<aug_t>*, ParallelTopologyCluster<aug_t>*> clusters;
    std::multimap<ParallelTopologyCluster<aug_t>*, ParallelTopologyCluster<aug_t>*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<ParallelTopologyCluster<aug_t>*, vertex_t> vertex_map;
    for (int i = 0; i < n; i++) vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < n; i++) clusters.insert({leaves[i].parent, &leaves[i]});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
        for (auto neighbor : leaf->neighbors) if (neighbor) std::cout << vertex_map[neighbor] << " ";
        std::cout << std::endl;
        bool in_map = false;
        for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
        if (parent && !in_map) next_clusters.insert({parent->parent, parent});
    }
    clusters.swap(next_clusters);
    next_clusters.clear();
    while (!clusters.empty()) {
        std::cout << "======================= NEXT LEVEL =======================" << std::endl;
        for (auto entry : clusters) {
            auto cluster = entry.second;
            auto parent = entry.first;
            std::cout << "Cluster: " << cluster << " Parent: " << parent << std::endl;
            bool in_map = false;
            for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
            if (parent && !in_map) next_clusters.insert({parent->parent, parent});
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
}

TEST(ParallelTopologyTreeSuite, batch_incremental_linkedlist_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
    }
}

TEST(ParallelTopologyTreeSuite, batch_incremental_binarytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
    }
}

TEST(ParallelTopologyTreeSuite, batch_incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::random_degree3_benchmark(n, rand());
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) return;
            tree.batch_link(batch.edges);
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after batch of links.";
        }
    }
}

TEST(ParallelTopologyTreeSuite, batch_decremental_linkedlist_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

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

TEST(ParallelTopologyTreeSuite, batch_decremental_binarytree_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

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

TEST(ParallelTopologyTreeSuite, batch_decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        vertex_t k = 16;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        ParallelTopologyTree<int> tree(n, k, qt, f, 0, 0);

        auto update_sequence = dynamic_tree_benchmark::random_degree3_benchmark(n, rand());
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
