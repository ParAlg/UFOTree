#include <gtest/gtest.h>
#include <unordered_set>
#include "../include/ufo_tree.h"


template<typename v_t, typename e_t>
bool UFOTree<v_t, e_t>::is_valid() {
    std::unordered_map<Cluster*,int> clusters;
    std::unordered_map<Cluster*,int> next_clusters;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert({&this->leaves[i], 0});
    for (auto entry : clusters) { // Ensure that every pair of incident vertices are in the same component
        bool neighbors_connected = true;
        auto leaf = entry.first;
        FOR_ALL_NEIGHBORS(leaf, [&](Cluster* neighbor, e_t _) {
            if (leaf->get_root() != neighbor->get_root()) neighbors_connected = false;
        });
        if (!neighbors_connected) return false;
    }
    while (!clusters.empty()) {
        for (auto entry : clusters) {
            auto cluster = entry.first;
            if (cluster->fanout != entry.second) return false; // Ensure the fanout field is correct
            bool children_good = true;
            FOR_ALL_CHILDREN(cluster, [&](Cluster* child) { // Ensure all children point back
                if (child->parent != cluster) children_good = false;
            });
            if (!children_good) return false;
            bool neighbors_good = true;
            FOR_ALL_NEIGHBORS(cluster, [&](Cluster* neighbor, e_t _) { // Ensure all neighbors point back
                if (!neighbor->contains_neighbor(cluster)) neighbors_good = false;
            });
            if (!neighbors_good) return false;
            if (cluster->degree <= 3 && !cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->degree == 1) {
                    if (cluster->neighbors[0]->degree > 2) return false;
                    else if (!cluster->neighbors[0]->contracts()) return false;
                } else if (cluster->degree == 2) {
                    for (auto neighborp : cluster->neighbors) {
                        auto neighbor = UNTAG(neighborp);
                        if (neighbor && neighbor->degree < 3 && !neighbor->contracts()) return false;
                    }
                } else if (cluster->degree >= 3) {
                    bool neighbors_deg_1 = false;
                    FOR_ALL_NEIGHBORS(cluster, [&](Cluster* neighbor, e_t _) {
                        if (neighbor->degree == 1) neighbors_deg_1 = true;
                    });
                    if (neighbors_deg_1) return false;
                }
            }
            if (cluster->parent) { // Get next level
                if (next_clusters.find(cluster->parent) != next_clusters.end())
                    next_clusters[cluster->parent]++;
                else
                    next_clusters.insert({cluster->parent, 1});
            }
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

template<typename v_t, typename e_t>
void UFOTree<v_t, e_t>::print_tree() {
    std::multimap<Cluster*, Cluster*> clusters;
    std::multimap<Cluster*, Cluster*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<Cluster*, vertex_t> vertex_map;
    for (int i = 0; i < this->leaves.size(); i++) vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert({leaves[i].parent, &leaves[i]});
    for (auto entry : clusters) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << vertex_map[leaf] << "\t " << leaf << " Parent " << parent << " Neighbors: ";
        FOR_ALL_NEIGHBORS(leaf, [&](Cluster* neighbor, e_t _) {
            std::cout << vertex_map[UNTAG(neighbor)] << " ";
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
            std::cout << "Cluster: " << cluster << " Parent: " << parent << std::endl;
            bool in_map = false;
            for (auto entry : next_clusters) if (entry.second == parent) in_map = true;
            if (parent && !in_map) next_clusters.insert({parent->parent, parent});
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
}

TEST(UFOTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(UFOTreeSuite, incremental_star_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << 0 << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, incremental_random_correctness_test) {
    int num_trials = 100;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        UFOTree<int, int> tree(n);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
            }
        }
    }
}

TEST(UFOTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_FALSE(tree.connected(i,i+1)) << "Vertex " << i << " and " << i+1 << " connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_FALSE(tree.connected(i,2*i+1)) << "Vertex " << i << " and " << 2*i+1 << " connected.";
        ASSERT_FALSE(tree.connected(i,2*i+2)) << "Vertex " << i << " and " << 2*i+2 << " connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_FALSE(tree.connected((n-1)/2,n-1)) << "Vertex " << (n-1)/2 << " and " << n-1 << " connected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(UFOTreeSuite, decremental_star_correctness_test) {
    vertex_t n = 256;
    UFOTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(0,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(0,i+1);
        ASSERT_FALSE(tree.connected(0,i+1)) << "Vertex " << 0 << " and " << i+1 << " connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << 0 << " and " << i+1 << ".";
    }
}

TEST(UFOTreeSuite, decremental_random_correctness_test) {
    int num_trials = 100;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        UFOTree<int, int> tree(n);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges[links++] = {u,v};
            }
        }
        for (auto edge : edges) {
            auto u = edge.first;
            auto v = edge.second;
            tree.cut(u,v);
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << u << " and " << v << ".";
        }
    }
}

// Query Tests


TEST(UFOTreeQuerySuite, LinkedListQueryTest){
  std::vector<int> test_vals = {10, 100, 1000};
  srand(time(NULL));
  int seed = 1; 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      QueryType q = PATH;
      auto f = [](int x, int y){return std::min(x,y);};
      int im = std::numeric_limits<int>::max();
      UFOTree<int, int> tree(n, q, f, f, im, im, im, im);
      vertex_t u = rand() % (n-1), v = rand() % n; if(v == u){v++;}
      //vertex_t u = 0, v = 7;
      if(v < u) std::swap(u,v);
      /*std::cout << "Seed: " << seed << "\n";
    std::cout << "u = " << u << "v = " << v << "\n";*/
      std::vector<std::tuple<int,int,int>> edges;
      int min_edge_val = std::numeric_limits<int>::max();
      for(int i = 0; i < n-1; i++){
        int new_edge = rand() % 100;
        tree.link(i, i+1, new_edge);
        if(i >= u && i < v) min_edge_val = std::min(min_edge_val, new_edge);
        edges.push_back({i, i+1, new_edge});
      }

      // Test return of min_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != min_edge_val){
        tree.print_tree();
        std::cout << "================ EDGES ================\n";
        for(auto edge : edges){
          std::cout << std::get<0>(edge) << " " << std::get<1>(edge) << " " << std::get<2>(edge) << "\n";
        }
        std::cout << "=======================================\n";
        std::cout << "Seed: " << seed << "\n";
        std::cout << "u = " << u << ", v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, min_edge_val);
    }
  }
}

TEST(UFOTreeQuerySuite, BinaryTreeQueryTest){
  std::vector<int> test_vals = {7, 31, 127, 1023};
  srand(time(NULL));
  int seed = 1; 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      QueryType q = PATH;
      auto f = [](int x, int y){return std::min(x,y);};
      int im = std::numeric_limits<int>::max();
      UFOTree<int, int> tree(n, q, f, f, im, im, im, im);
      std::vector<int> path; 
      for(int i = 0; i < n - 1; i = (2*i) + 2){
        path.push_back(i);
      }

      vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();if(lower == upper) lower++;
      if(lower < upper) std::swap(lower, upper);

      auto u = path[upper], v = path[lower];
      //vertex_t u = 4, v = 7;
      //if(v < u) std::swap(u,v);
      /*std::cout << "Seed: " << seed << "\n";
      std::cout << "u = " << u << "v = " << v << "\n";*/ 
      int j = 0;
      int min_edge_val = std::numeric_limits<int>::max();
      for(int i = 0; i < (n/2); i++){ 
        auto new_edge = rand() % 100, new_edge2 = rand() % 100; 
        tree.link(i, (2*i) + 1, new_edge);
        tree.link(i, (2*i) + 2, new_edge2);
        if(i == path[j]){
          if(i >= u && i < v){
            min_edge_val = std::min(min_edge_val, new_edge2);
          }
          j++;
        }
      }

      // Test return of min_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != min_edge_val){
        tree.print_tree();
        std::cout << seed << "\n";
        std::cout << "u = " << u << "v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, min_edge_val);
    }
  }
}

TEST(UFOTreeQuerySuite, HighDegreeQuery){
  // Possible experiment setup:
  // 1) Star graph with a superunary node that continues until the root of the tree.
  // 2) Query with one node that contracts early and one that stays until the end.
}
