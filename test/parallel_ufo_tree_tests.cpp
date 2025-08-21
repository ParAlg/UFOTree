#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "parallel_ufo_tree.h"
#include "../benchmark/dynamic_trees/benchmark.h"
#include "../benchmark/parallel/parallel_benchmark.h"


using namespace dgbs;

template <typename aug_t>
bool ParallelUFOTree<aug_t>::is_valid(parlay::sequence<std::pair<int, int>>& edges) {
    // Ensure that the level 0 forest matches the input edge list
    auto dir_edges = parlay::tabulate(2*edges.size(), [&] (size_t i) {
        if (i % 2 == 0) return std::make_pair(edges[i/2].first, edges[i/2].second);
        return std::make_pair(edges[i/2].second, edges[i/2].first);
    });
    std::atomic<bool> all_edges_in_tree = true;
    parlay::parallel_for(0, dir_edges.size(), [&] (size_t i) {
        if (!leaves[dir_edges[i].first].contains_neighbor(&leaves[dir_edges[i].second]))
            all_edges_in_tree = false;
    });
    if (!all_edges_in_tree) {
        std::cerr << "NOT ALL EDGES ARE PRESENT IN THE TREE." << std::endl;
        return false;
    }
    std::set<std::pair<int, int>> edge_set;
    for (auto edge : dir_edges) edge_set.insert(edge);
    std::atomic<bool> no_extra_edges_in_tree = true;
    parlay::parallel_for(0, leaves.size(), [&] (size_t i) {
        leaves[i].for_all_neighbors([&] (auto neighbor) {
            int index = neighbor - &leaves[0];
            if (!edge_set.contains(std::make_pair(i, index)))
                no_extra_edges_in_tree = false;
        });
    });
    if (!no_extra_edges_in_tree) {
        std::cerr << "THERE ARE EXTRA EDGES PRESENT IN THE TREE." << std::endl;
        return false;
    }
    // Validate the tree based on its level 0 contents
    std::unordered_set<Cluster*> clusters;
    std::unordered_set<Cluster*> next_clusters;
    for (int i = 0; i < leaves.size(); i++) { // Ensure that every pair of incident vertices are in the same component
        std::atomic<bool> connectivity_correct = true;
        leaves[i].for_all_neighbors([&] (auto neighbor) {
            if (leaves[i].get_root() != neighbor->get_root()) connectivity_correct = false;
        });
        if (!connectivity_correct) {
            std::cerr << "CONNECTIVITY INCORRECT" << std::endl;
            return false;
        }
    }
    for (int i = 0; i < leaves.size(); i++) clusters.insert(&leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            std::atomic<bool> pointers_bidirectional = true; // Ensure all neighbors also point back
            cluster->for_all_neighbors([&] (auto neighbor) {
                if (!neighbor->contains_neighbor(cluster)) pointers_bidirectional = false;
            });
            if (!pointers_bidirectional) {
                std::cerr << "NEIGHBOR DOESN'T POINT BACK, CLUSTER " << cluster << std::endl;
                return false;
            }
            if (cluster->get_degree() <= 3 && !cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    if (cluster->get_neighbor()->get_degree() > 2 || !cluster->get_neighbor()->contracts()) {
                        std::cerr << "NON-MAXIMAL CONTRACTION, DEG 1 CLUSTER " << cluster << std::endl;
                        return false;
                    }
                } else if (cluster->get_degree() == 2) {
                    std::atomic<bool> deg2_maximal = true;
                    cluster->for_all_neighbors([&] (auto neighbor) {
                        if (neighbor->get_degree() < 3 && !neighbor->contracts()) deg2_maximal = false;
                    });
                    if (!deg2_maximal) {
                        std::cerr << "NON-MAXIMAL CONTRACTION, DEG 2 CLUSTER " << cluster << std::endl;
                        return false;
                    }
                } else if (cluster->get_degree() >= 3) {
                    std::atomic<bool> deg3_maximal = true;
                    cluster->for_all_neighbors([&] (auto neighbor) {
                        if (neighbor->get_degree() == 1 && neighbor->parent != cluster->parent) deg3_maximal = false;
                    });
                    if (!deg3_maximal) {
                        std::cerr << "NON-MAXIMAL CONTRACTION, DEG 3+ CLUSTER " << cluster << std::endl;
                        return false;
                    }
                }
            }
            for (auto cluster : clusters) { // Ensure no partner fields are still set
                if (cluster->partner) {
                    std::cerr << "PARTNER FIELD ERRONEOUSLY SET, CLUSTER " << cluster << std::endl;
                    return false;
                }
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
        // leaf->print_neighbors();
        ParallelUFOCluster<aug_t>::ufo_pam_set::foreach_seq(leaf->neighbors, [&] (auto neighbor) {
            std::cout << neighbor - &leaves[0] << " ";
        });
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
    ufo_pam_set::foreach_seq(neighbors, [&] (auto neighbor) { std::cout << neighbor << " "; });
    std::cout << std::endl;
}


parlay::internal::timer timer0("");
parlay::internal::timer timer1("");
parlay::internal::timer timer2("");
parlay::internal::timer timer3("");
parlay::internal::timer timer4("");
parlay::internal::timer timer5("");

parlay::internal::timer subtimer1("");
parlay::internal::timer subtimer2("");
parlay::internal::timer subtimer3("");
parlay::internal::timer subtimer4("");

extern int command_line_n;
extern int command_line_k;
extern int command_line_num_trials;
extern int command_line_seed;

TEST(ParallelUFOTreeSuite, batch_incremental_linkedlist_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            edges = parlay::append(edges, batch.edges);
            if (!tree.is_valid(edges)) {
                tree.print_tree();
                FAIL() << "Tree invalid after batch of links.";
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_star_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            edges = parlay::append(edges, batch.edges);
            if (!tree.is_valid(edges)) {
                tree.print_tree();
                FAIL() << "Tree invalid after batch of links.";
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_binarytree_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            edges = parlay::append(edges, batch.edges);
            if (!tree.is_valid(edges)) {
                tree.print_tree();
                FAIL() << "Tree invalid after batch of links.";
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_karytree_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;
    vertex_t fanout = 8;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n, seed, fanout);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            edges = parlay::append(edges, batch.edges);
            if (!tree.is_valid(edges)) {
                tree.print_tree();
                FAIL() << "Tree invalid after batch of links.";
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_incremental_random_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type != INSERT) break;
            tree.batch_link(batch.edges);
            edges = parlay::append(edges, batch.edges);
            if (!tree.is_valid(edges)) {
                tree.print_tree();
                FAIL() << "Tree invalid after batch of links.";
            }
        }
    }
}

// ===================================================
// ================ DECREMENTAL TESTS ================
// ===================================================

TEST(ParallelUFOTreeSuite, batch_decremental_linkedlist_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
                edges = parlay::append(edges, batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                for (auto edge : batch.edges)
                    edges = parlay::remove(edges, edge);
                if (!tree.is_valid(edges)) {
                    tree.print_tree();
                    FAIL() << "Tree invalid after batch of links.";
                }
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_star_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
                edges = parlay::append(edges, batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                for (auto edge : batch.edges)
                    edges = parlay::remove(edges, edge);
                if (!tree.is_valid(edges)) {
                    tree.print_tree();
                    FAIL() << "Tree invalid after batch of links.";
                }
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_binarytree_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
                edges = parlay::append(edges, batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                for (auto edge : batch.edges)
                    edges = parlay::remove(edges, edge);
                if (!tree.is_valid(edges)) {
                    tree.print_tree();
                    FAIL() << "Tree invalid after batch of links.";
                }
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_karytree_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;
    vertex_t fanout = 8;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n, seed, fanout);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
                edges = parlay::append(edges, batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                for (auto edge : batch.edges)
                    edges = parlay::remove(edges, edge);
                if (!tree.is_valid(edges)) {
                    tree.print_tree();
                    FAIL() << "Tree invalid after batch of links.";
                }
            }
        }
    }
}

TEST(ParallelUFOTreeSuite, batch_decremental_random_correctness_test) {
    vertex_t n = command_line_n > 0 ? command_line_n : 256;
    vertex_t k = command_line_k > 0 ? command_line_k : 16;
    int num_trials = command_line_num_trials > 0 ? command_line_num_trials : 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        ParallelUFOTree<> tree(n, k);
        long seed = command_line_seed != -1 ? command_line_seed : seeds[trial];
        std::cout << "SEED: " << seed << std::endl;

        auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n, seed);
        auto batches = parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k);
        parlay::sequence<std::pair<int, int>> edges;

        for (auto batch : batches) {
            if (batch.type == INSERT) {
                tree.batch_link(batch.edges);
                edges = parlay::append(edges, batch.edges);
            } else {
                tree.batch_cut(batch.edges);
                for (auto edge : batch.edges)
                    edges = parlay::remove(edges, edge);
                if (!tree.is_valid(edges)) {
                    tree.print_tree();
                    FAIL() << "Tree invalid after batch of links.";
                }
            }
        }
    }
}
