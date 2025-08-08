#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "parallel_ufo_tree.h"
#include "../benchmark/dynamic_trees/benchmark.h"
#include "../benchmark/parallel/parallel_benchmark.h"


using namespace dgbs;

template <typename aug_t>
bool ParallelUFOTree<aug_t>::is_valid() {
    std::unordered_set<Cluster*> clusters;
    std::unordered_set<Cluster*> next_clusters;
    for (int i = 0; i < leaves.size(); i++) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaves[i].neighbors) // This ensures all connectivity is correct by transitivity
            if (leaves[i].get_root() != neighbor->get_root()) {
                std::cerr << "CONNECTIVITY INCORRECT" << std::endl;
                return false;
            }
    for (int i = 0; i < leaves.size(); i++) clusters.insert(&leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (!neighbor->contains_neighbor(cluster)) {
                    std::cerr << "NEIGHBOR DOESN'T POINT BACK, CLUSTER " << cluster << std::endl;
                    return false;
                }
            if (cluster->get_degree() <= 3 && !cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    if (cluster->get_neighbor()->get_degree() > 2) {
                        std::cerr << "NON-MAXIMAL CONTRACTION, DEG 1 CLUSTER " << cluster << std::endl;
                        return false;
                    }
                    else if (!cluster->get_neighbor()->contracts()) {
                        std::cerr << "NON-MAXIMAL CONTRACTION, DEG 1 CLUSTER " << cluster << std::endl;
                        return false;
                    }
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor->get_degree() < 3 && !neighbor->contracts()) {
                            std::cerr << "NON-MAXIMAL CONTRACTION, DEG 2 CLUSTER " << cluster << std::endl;
                            return false;
                        }
                } else if (cluster->get_degree() >= 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && neighbor->get_degree() < 2) {
                            std::cerr << "NON-MAXIMAL CONTRACTION, DEG 3+ CLUSTER " << cluster << std::endl;
                            return false;
                        }
                }
            }
            for (auto cluster : clusters) // Ensure no partner fields are still set
                if (cluster->partner) {
                    std::cerr << "PARTNER FIELD ERRONEOUSLY SET, CLUSTER " << cluster << std::endl;
                    return false;
                }
            if (cluster->parent) // Get next level
                next_clusters.insert(cluster->parent);
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

template <typename aug_t>
int ParallelUFOTree<aug_t>::get_height(vertex_t v) {
    // return forests.size();
    return 0;
}

template <typename aug_t>
void ParallelUFOTree<aug_t>::print_tree() {
    std::multimap<Cluster*, Cluster*> clusters;
    std::multimap<Cluster*, Cluster*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<Cluster*, vertex_t> vertex_map;
    for (int i = 0; i < leaves.size(); i++) vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < leaves.size(); i++) clusters.insert({leaves[i].parent, &leaves[i]});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
        for (auto neighbor : leaf->neighbors) std::cout << vertex_map[neighbor] << " ";
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
            std::cout << "Cluster: " << cluster << " Parent: " << parent << " Neighbors: ";
            cluster->print_neighbors();
            bool in_map = false;
            for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
            if (parent && !in_map) next_clusters.insert({parent->parent, parent});
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
}

template <typename aug_t>
void ParallelUFOCluster<aug_t>::print_neighbors() {
    for (auto neighbor : neighbors) std::cout << neighbor << " ";
    std::cout << std::endl;
}

TEST(ParallelUFOTreeSuite, batch_incremental_linkedlist_correctness_test) {
    int num_trials = 100;
    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 6;
        vertex_t k = 1;
        ParallelUFOTree<> tree(n, k);
        long seed = seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            // tree.print_tree();
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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
        ParallelUFOTree<> tree(n, k);

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
