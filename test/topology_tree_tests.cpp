#include <gtest/gtest.h>
#include <limits>
#include <unordered_set>
#include "../include/topology_tree.h"


template<typename v_t, typename e_t>
bool TopologyTree<v_t, e_t>::is_valid() {
    std::unordered_set<TopologyCluster<v_t, e_t>*> clusters;
    std::unordered_set<TopologyCluster<v_t, e_t>*> next_clusters;
    for (auto leaf : leaves) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor : leaf.neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor && leaf.get_root() != neighbor->get_root()) return false;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(&this->leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor && !neighbor->contains_neighbor(cluster)) return false;
            if (!cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts()) return false;
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 3) return false;
                } else if (cluster->get_degree() == 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 2) return false;
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

template<typename v_t, typename e_t>
void TopologyTree<v_t, e_t>::print_tree() {
    std::multimap<TopologyCluster<v_t, e_t>*, TopologyCluster<v_t, e_t>*> clusters;
    std::multimap<TopologyCluster<v_t, e_t>*, TopologyCluster<v_t, e_t>*> next_clusters;
    std::cout << "========================= LEAVES =========================" << std::endl;
    std::unordered_map<TopologyCluster<v_t, e_t>*, vertex_t> vertex_map;
    for (int i = 0; i < this->leaves.size(); i++) vertex_map.insert({&leaves[i], i});
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert({leaves[i].parent, &leaves[i]});
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

TEST(TopologyTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    TopologyTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    TopologyTree<int, int> tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(TopologyTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; ++trial) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; ++trial) {
        vertex_t n = 256;
        TopologyTree<int, int> tree(n);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        std::vector<int> vertex_degrees(n,0);
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (vertex_degrees[u] >= 3) continue;
            if (vertex_degrees[v] >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
                vertex_degrees[u]++;
                vertex_degrees[v]++;
            }
        }
    }
}

TEST(TopologyTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    TopologyTree<int, int> tree(n);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    TopologyTree<int, int> tree(n);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(TopologyTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; ++trial) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; ++trial) {
        vertex_t n = 256;
        TopologyTree<int, int> tree(n);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        std::vector<int> vertex_degrees(n,0);
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (vertex_degrees[u] >= 3) continue;
            if (vertex_degrees[v] >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges[links++] = {u,v};
                vertex_degrees[u]++;
                vertex_degrees[v]++;
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

TEST(TopologyTreeSuite, subtree_query_test) {
    vertex_t n = 256;
    QueryType qt = SUBTREE;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 1;
    TopologyTree<int, int> tree(n, qt, f, f, id, id, d, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t v = 0; v <= i+1; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i+2-v);
            ASSERT_EQ(tree.subtree_query(v, (v <= i ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
    for (vertex_t i = n-1; i > 0; i--) {
        tree.cut(i-1,i);
        for (vertex_t v = 0; v < i; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i-v);
            ASSERT_EQ(tree.subtree_query(v, (v < i-1 ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
}

TEST(TopologyTreeSuite, path_query_test) {
    vertex_t n = 64;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 0;
    TopologyTree<int, int> tree(n, qt, f, f, id, id, d, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1,1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_EQ(tree.path_query(u,v), v-u);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_EQ(tree.path_query(u,v), v-u);
    }
}

TEST(TopologyTreeQuerySuite, LinkedListQueryTest){
  std::vector<int> test_vals = {10, 100, 1000};
  int num_trials = 1;
  int seeds[num_trials];
  srand(time(NULL));
  for (int trial = 0; trial < num_trials; ++trial) seeds[trial] = rand();
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      int seed = seeds[trial];
      srand(seed);
      QueryType qt = PATH;
      auto f = [](int x, int y)->int {return std::min(x,y);};
      int id = std::numeric_limits<int>::max();
      int d = std::numeric_limits<int>::max();
      TopologyTree<int, int> tree(n, qt, f, f, id, id, d, d);
      vertex_t u = rand() % (n-1), v = rand() % n; if(v == u){v++;}
      if(v < u) std::swap(u,v);
      int min_edge_val = std::numeric_limits<int>::max();
      for(int i = 0; i < n-1; i++){
        int new_edge = rand() % 100;
        tree.link(i, i+1, new_edge);
        if(i >= u && i < v) min_edge_val = std::min(min_edge_val, new_edge);
      }

      // Test return of min_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != min_edge_val){
        tree.print_tree();
        std::cout << "SEED: " << seed << std::endl;
        std::cout << "u = " << u << ", v = " << v << std::endl;
      }
      ASSERT_EQ(returned_query, min_edge_val);
    }
  }
}

TEST(TopologyTreeQuerySuite, BinaryTreeQueryTest){
  std::vector<int> test_vals = {15, 127, 1023};
  int num_trials = 1;
  int seeds[num_trials];
  srand(time(NULL));
  for (int trial = 0; trial < num_trials; ++trial) seeds[trial] = rand();
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      int seed = seeds[trial];
      srand(seed);
      QueryType qt = PATH;
      auto f = [](int x, int y)->int {return std::min(x,y);};
      int id = std::numeric_limits<int>::max();
      int d = std::numeric_limits<int>::max();
      TopologyTree<int, int> tree(n, qt, f, f, id, id, d, d);
      std::vector<int> path; 
      for(int i = 0; i < n - 1; i = (2*i) + 2){
        path.push_back(i);
      }
      vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();if(lower == upper) lower++;
      if(lower < upper) std::swap(lower, upper);
      auto u = path[upper], v = path[lower];
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
        std::cout << "SEED: " << seed << std::endl;
        std::cout << "u = " << u << ", v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, min_edge_val);
    }
  }
}
