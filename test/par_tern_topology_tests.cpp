#pragma once
#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "util.h"
#include "types.h"
#include "parallel_topology_tree_ternarized.h"
#include "../benchmark/parallel/parallel_benchmark.h"
#include "../benchmark/dynamic_trees/benchmark.h"
  

namespace ufo {
  static std::vector<UpdateBatchWithWeights> convert_to_weighted_update_batch(std::vector<UpdateBatch> update_batches){
    // Generate weighted Update Sequence for RC trees (only in the case of insertion)
      std::vector<UpdateBatchWithWeights> weighted_update_sequence;
      weighted_update_sequence.resize(update_batches.size());
      // Initialzing same update sequence with default val edge weights
      for (int i = 0; i <  update_batches.size(); i++) {
        auto updateBatch = update_batches[i];
        auto edges = updateBatch.edges;
        weighted_update_sequence[i].type = updateBatch.type;
        if (updateBatch.type == INSERT) {
          for (int k = 0; k < edges.size(); k++) {
            weighted_update_sequence[i].insert_edges.push_back(std::make_tuple(edges[k].first, edges[k].second, 0));
          }
        }
        else {
          weighted_update_sequence[i].delete_edges = edges;
        }
      }
      return weighted_update_sequence;
  }
  template <typename aug_t>
  bool ParallelTopologyTreeTernarized<aug_t>::is_valid() {
    std::unordered_set<ParallelTopologyCluster<aug_t>*> clusters;
    std::unordered_set<ParallelTopologyCluster<aug_t>*> next_clusters;
    for (int i = 0; i < n; i++)                 // Ensure that every pair of incident vertices are in the same component
      for (auto neighbor : leaves[i].neighbors) // This ensures all connectivity is correct by transitivity
        if (neighbor && leaves[i].get_root() != neighbor->get_root())
        {
          std::cerr << "INCORRECT CONNECTIVITY." << std::endl;
          return false;
        }
    for (int i = 0; i < n; i++)
      clusters.insert(&leaves[i]);
    while (!clusters.empty())
    {
      for (auto cluster : clusters)
      {
        for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
          if (neighbor && !neighbor->contains_neighbor(cluster))
          {
            std::cerr << "INCORRECT ADJACENCY." << std::endl;
            return false;
          }
        if (!cluster->contracts())
        { // Ensure maximality of contraction
          if (cluster->get_degree() == 1)
          {
            for (auto neighbor : cluster->neighbors)
              if (neighbor && !neighbor->contracts())
              {
                std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                return false;
              }
          }
          else if (cluster->get_degree() == 2)
          {
            for (auto neighbor : cluster->neighbors)
              if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 3)
              {
                std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                return false;
              }
          }
          else if (cluster->get_degree() == 3)
          {
            for (auto neighbor : cluster->neighbors)
              if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 2)
              {
                std::cerr << "CONTRACTIONS NOT MAXIMAL." << std::endl;
                return false;
              }
          }
        }
        if (cluster->parent)
          next_clusters.insert(cluster->parent); // Get next level
      }
      clusters.swap(next_clusters);
      next_clusters.clear();
    }
    return true;
  }

  template <typename aug_t>
  int ParallelTopologyTreeTernarized<aug_t>::get_height(vertex_t v)
  {
    int height = 0;
    ParallelTopologyCluster<aug_t> *curr = &leaves[v];
    while (curr)
    {
      height++;
      curr = curr->parent;
    }
    return height;
  }

  template <typename aug_t>
  void ParallelTopologyTreeTernarized<aug_t>::print_tree()
  {
    std::multimap<ParallelTopologyCluster<aug_t> *, ParallelTopologyCluster<aug_t> *> clusters;
    std::multimap<ParallelTopologyCluster<aug_t> *, ParallelTopologyCluster<aug_t> *> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<ParallelTopologyCluster<aug_t> *, vertex_t> vertex_map;
    for (int i = 0; i < n; i++)
      vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < n; i++)
      clusters.insert({leaves[i].parent, &leaves[i]});
    for (auto entry : clusters)
    {
      auto leaf = entry.second;
      auto parent = entry.first;
      std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
      for (auto neighbor : leaf->neighbors)
        if (neighbor)
          std::cout << vertex_map[neighbor] << " ";
      std::cout << std::endl;
      bool in_map = false;
      for (auto entry : next_clusters)
        if (entry.second == parent)
          in_map = true;
      if (parent && !in_map)
        next_clusters.insert({parent->parent, parent});
    }
    clusters.swap(next_clusters);
    next_clusters.clear();
    while (!clusters.empty())
    {
      std::cout << "======================= NEXT LEVEL =======================" << std::endl;
      for (auto entry : clusters)
      {
        auto cluster = entry.second;
        auto parent = entry.first;
        std::cout << "Cluster: " << cluster << " Parent: " << parent << std::endl;
        bool in_map = false;
        for (auto entry : next_clusters)
          if (entry.second == parent)
            in_map = true;
        if (parent && !in_map)
          next_clusters.insert({parent->parent, parent});
      }
      clusters.swap(next_clusters);
      next_clusters.clear();
    }
  }
  TEST(ParallelTernarizedTopologySuite, test_constructor)
  {
    ParallelTopologyTreeTernarized<int> t(3, 3);
    // t.print_tree();
  } 

  TEST(ParallelTernarizedTopologySuite, batch_incremental_linkedlist_correctness_test)
  {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials =  1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++)
    {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches)
      {
        if (batch.type != INSERT)
          break;
        tree.batch_link(batch.insert_edges);
        if (!tree.is_valid())
        {
          tree.print_tree();
          FAIL() << "Tree invalid after batch of links.";
        }
      }
    }
  }
  
  TEST(ParallelTernarizedTopologySuite, batch_incremental_star_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 1;
    int num_trials =  1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));
      
      for (auto batch : batches)
      {
        if (batch.type != INSERT)
          break;
          tree.batch_link(batch.insert_edges);
        if (!tree.is_valid()) {
          tree.print_tree();
          FAIL() << "Tree invalid after batch of links.";
        }
      }
    }
  }

  TEST(ParallelTernarizedTopologySuite, batch_incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches){
        if (batch.type != INSERT)
          break;
        tree.batch_link(batch.insert_edges);
        if (!tree.is_valid()){
          tree.print_tree();
          FAIL() << "Tree invalid after batch of links.";
        }
      }
    }
  }

  TEST(ParallelTernarizedTopologySuite, batch_incremental_karytree_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;
    vertex_t fanout = 8;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n, seed, fanout);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type != INSERT)
          break;
        tree.batch_link(batch.insert_edges);
        if (!tree.is_valid()) {
          tree.print_tree();
          FAIL() << "Tree invalid after batch of links.";
        }
      }
    }
  }

  TEST(ParallelTernarizedTopologySuite, batch_incremental_random_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));
      parlay::sequence<std::pair<int, int>> edges;

      for (auto batch : batches) {
        if (batch.type != INSERT)
          break;
        tree.batch_link(batch.insert_edges);
        if (!tree.is_valid()) {
          tree.print_tree();
          FAIL() << "Tree invalid after batch of links.";
        }
      }
    }
  }
  // ===================================================
  // ================ DECREMENTAL TESTS ================
  // ===================================================

  TEST(ParallelTernarizedTopologySuite, batch_decremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::linked_list_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type == INSERT) {
          tree.batch_link(batch.insert_edges);
        }
        else {
          tree.batch_cut(batch.delete_edges);
          
          if (!tree.is_valid()) {
            tree.print_tree();
            FAIL() << "Tree invalid after batch of links.";
          }
        }
      }
    }
  }
  
  TEST(ParallelTernarizedTopologySuite, batch_decremental_star_correctness_test) { 
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::star_benchmark(n, seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type == INSERT) {
          tree.batch_link(batch.insert_edges);
        }
        else {
          tree.batch_cut(batch.delete_edges);
          
          if (!tree.is_valid()) {
            tree.print_tree();
            FAIL() << "Tree invalid after batch of links.";
          }
        }
      }
    }
  }
  TEST(ParallelTernarizedTopologySuite, batch_decremental_binarytree_correctness_test) {

    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::binary_tree_benchmark(n,seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type == INSERT) {
          tree.batch_link(batch.insert_edges);
        }
        else {
          tree.batch_cut(batch.delete_edges);
          
          if (!tree.is_valid()) {
            tree.print_tree();
            FAIL() << "Tree invalid after batch of links.";
          }
        }
      }
    } 
}

  TEST(ParallelTernarizedTopologySuite, batch_decremental_karytree_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    vertex_t fanout = 8;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::k_ary_tree_benchmark_helper(n,seed,fanout);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type == INSERT) {
          tree.batch_link(batch.insert_edges);
        }
        else {
          tree.batch_cut(batch.delete_edges);
          
          if (!tree.is_valid()) {
            tree.print_tree();
            FAIL() << "Tree invalid after batch of links.";
          }
        }
      }
    } 
  }

  TEST(ParallelTernarizedTopologySuite, batch_decremental_random_correctness_test) {
    vertex_t n = 256;
    vertex_t k = 16;
    int num_trials = 1;

    long seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++)
      seeds[trial] = rand();

    for (int trial = 0; trial < num_trials; trial++) {
      ParallelTopologyTreeTernarized<int> tree(n, k);
      long seed = seeds[trial];
      std::cout << "SEED: " << seed << std::endl;

      auto update_sequence = dynamic_tree_benchmark::random_unbounded_benchmark(n,seed);
      auto batches = convert_to_weighted_update_batch(parallel_dynamic_tree_benchmark::convert_updates_to_batches(update_sequence, k));

      for (auto batch : batches) {
        if (batch.type == INSERT) {
          tree.batch_link(batch.insert_edges);
        }
        else {
          tree.batch_cut(batch.delete_edges);
          
          if (!tree.is_valid()) {
            tree.print_tree();
            FAIL() << "Tree invalid after batch of links.";
          }
        }
      }
    }
  }
}