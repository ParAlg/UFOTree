#include <gtest/gtest.h>
#include <unordered_set>
#include <parlay/sequence.h>
#include <parlay/primitives.h>
#include "parallel_ufo_tree.h"
#include "../benchmark/dynamic_trees/benchmark.h"
#include "../benchmark/parallel/parallel_benchmark.h"


using namespace dgbs;
using Cluster = ParallelUFOCluster<empty_t>;
using UFOTree = ParallelUFOTree<empty_t>;
using allocator = parlay::type_allocator<Cluster>;

TEST(ParallelMISSuite, mis_test) {
    int num_trials = 1000000;
    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
        std::cout << "TRIAL " << trial << std::endl;
        int n = 100;
        long seed = seeds[trial];
        parlay::sequence<Cluster> clusters(n);

        // Randomly permute the clusters and connect them as a linked list
        parlay::sequence<int> permutation = parlay::random_permutation(n, seed);
        parlay::parallel_for(0, n, [&] (size_t i) {
            if (i < n-1) {
                clusters[permutation[i]].insert_neighbor(&clusters[permutation[i+1]]);
                clusters[permutation[i+1]].insert_neighbor(&clusters[permutation[i]]);
            }
        });

        // Initialize all clusters as root clusters
        parlay::sequence<Cluster*> root_clusters = parlay::tabulate(n, [&] (size_t i) {
            return &clusters[permutation[i]];
        });

        // Call the parallel MIS algorithm to recluster the root clusters
        UFOTree::recluster_root_clusters(root_clusters);
        UFOTree::create_new_parents(root_clusters);

        // Check that all clusters have a parent
        std::atomic<Cluster*> parentless_cluster = nullptr;
        parlay::parallel_for(0, n, [&] (size_t i) {
            if (!root_clusters[i]->parent)
                parentless_cluster = root_clusters[i];
        });
        ASSERT_EQ(parentless_cluster, nullptr) << "ROOT CLUSTER " << parentless_cluster << " DID NOT GET A PARENT." << std::endl;

        // Check that all contractions are maximal
        std::atomic<Cluster*> non_maximal_cluster = nullptr;
        parlay::parallel_for(0, n, [&] (size_t i) {
            Cluster* cluster = root_clusters[i];
            if (!cluster->contracts()) {
                if (cluster->get_degree() == 1) {
                    if (cluster->get_neighbor()->get_degree() > 2)
                        non_maximal_cluster = cluster;
                    else if (!cluster->get_neighbor()->contracts())
                        non_maximal_cluster = cluster;
                } else if (cluster->get_degree() == 2) {
                    cluster->for_all_neighbors([&] (auto neighbor) {
                        if (neighbor->get_degree() < 3 && !neighbor->contracts())
                            non_maximal_cluster = cluster;
                    });
                }
            }
        });
        ASSERT_EQ(non_maximal_cluster, nullptr) << "ROOT CLUSTER " << non_maximal_cluster << " DOES NOT CONTRACT MAXIMALLY." << std::endl;
    }
}
