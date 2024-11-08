#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/top_tree.h"

TEST(TopTreeSuite, constructor_test) {
    TopTree<int> t(3, PATH, [] (int a, int b){return std::min(a,b);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    ASSERT_EQ(t.t.num_vertices, 3);
}

TEST(TopTreeSuite, basic_link_test) {
  int n = 100;
  TopTree<int> t(n, PATH, [] (int a, int b){return std::min(a,b);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
  for (int i = 0; i < n-1; ++i) {
    t.link(i,i+1,i);
    for (int j = 0; j <= i; ++j) ASSERT_TRUE(t.connected(j,j+1));
    for (int j = i+1; j < n-1; ++j) ASSERT_FALSE(t.connected(j,j+1)); 
  }
}

TEST(TopTreeSuite, testInsertDeleteRakeLinkedList) {
  std::vector<vertex_t> to_test({10, 100});
  for (auto n : to_test) {
    TopTree<int> tree(n);
    for (int i = 0; i < n-1; ++i) {
      tree.link(i,i+1,i);
      for (int j = 0; j <= i; ++j) ASSERT_TRUE(tree.connected(j,j+1));
      for (int j = i+1; j < n-1; ++j) ASSERT_FALSE(tree.connected(j,j+1)); 
    }
    for (int i = 0; i < n-1; i++) {
      tree.cut(i,i+1);
      for (int j = 0; j <= i; ++j) ASSERT_FALSE(tree.connected(j,j+1));
      for (int j = i+1; j < n-1; ++j) ASSERT_TRUE(tree.connected(j,j+1)); 
    }
  }
}

TEST(Top_tree_suite, testLinkedListMaxQueries){
  std::vector<int> test_vals = {10, 100, 1000};
  srand(time(NULL));
  int seed = rand(); 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      TopTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::max(x,y);}, std::numeric_limits<int>::max(), 0);
      vertex_t u = rand() % (n-1), v = rand() % n; if(v == u){v++;}
      //vertex_t u = 4, v = 7;
      if(v < u) std::swap(u,v);
      /*std::cout << "Seed: " << seed << "\n";
    std::cout << "u = " << u << "v = " << v << "\n";*/
      int max_edge_val = std::numeric_limits<int>::min();
      for(int i = 0; i < n-1; i++){
        int new_edge = rand() % 100;
        tree.link(i, i+1, new_edge);
        if(i >= u && i < v) max_edge_val = std::max(max_edge_val, new_edge);
      }

      // Test return of max_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != max_edge_val){
        std::cout << seed << "\n";
        std::cout << "u = " << u << "v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, max_edge_val);
    }
  }
}

TEST(Top_tree_suite, BinaryTreeMaxQueryTest){
  std::vector<int> test_vals = {7, 31, 127, 1023};
  srand(time(NULL));
  int seed = rand(); 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      TopTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::max(x,y);}, std::numeric_limits<int>::max(), 0);
      std::vector<int> path;
      for(int i = 0; i < n - 1; i = (2*i) + 2){
        path.push_back(i);
      }
      
      vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();if(lower == upper) lower++;
      if(lower < upper) std::swap(lower, upper);

      auto u = path[upper], v = path[lower];
      //ASSERT_TRUE(u < v);
      //vertex_t u = 4, v = 7;
      /*std::cout << "Seed: " << seed << "\n";
      std::cout << "u = " << u << "v = " << v << "\n";*/
      int j = 0;
      int max_edge_val = std::numeric_limits<int>::min();
      for(int i = 0; i < (n/2); i++){ 
        auto new_edge = rand() % 100, new_edge2 = rand() % 100; 
        tree.link(i, (2*i) + 1, new_edge);
        tree.link(i, (2*i) + 2, new_edge2);
        if(i == path[j]){
          if(i >= u && i < v){
            max_edge_val = std::max(max_edge_val, new_edge2);
          }
          j++;
        }
      }

      // Test return of max_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != max_edge_val){
        std::cout << seed << "\n";
        std::cout << "u = " << u << "v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, max_edge_val);
    }
  }
}

TEST(Top_tree_suite, testLinkedListMaxQueries){
  std::vector<int> test_vals = {10, 100, 1000};
  srand(time(NULL));
  int seed = rand(); 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      TopTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::max(x,y);}, std::numeric_limits<int>::max(), 0);
      vertex_t u = rand() % (n-1), v = rand() % n; if(v == u){v++;}
      //vertex_t u = 4, v = 7;
      if(v < u) std::swap(u,v);
      /*std::cout << "Seed: " << seed << "\n";
    std::cout << "u = " << u << "v = " << v << "\n";*/
      int max_edge_val = std::numeric_limits<int>::min();
      for(int i = 0; i < n-1; i++){
        int new_edge = rand() % 100;
        tree.link(i, i+1, new_edge);
        if(i >= u && i < v) max_edge_val = std::max(max_edge_val, new_edge);
      }

      // Test return of max_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != max_edge_val){
        std::cout << seed << "\n";
        std::cout << "u = " << u << "v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, max_edge_val);
    }
  }
}

TEST(Top_tree_suite, BinaryTreeMaxQueryTest){
  std::vector<int> test_vals = {7, 31, 127, 1023};
  srand(time(NULL));
  int seed = rand(); 
  srand(seed);
  int num_trials = 1;
  for(int n : test_vals){
    for(int trial = 0; trial < num_trials; ++trial){
      TopTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::max(x,y);}, std::numeric_limits<int>::max(), 0);
      std::vector<int> path;
      for(int i = 0; i < n - 1; i = (2*i) + 2){
        path.push_back(i);
      }
      
      vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();if(lower == upper) lower++;
      if(lower < upper) std::swap(lower, upper);

      auto u = path[upper], v = path[lower];
      //ASSERT_TRUE(u < v);
      //vertex_t u = 4, v = 7;
      /*std::cout << "Seed: " << seed << "\n";
      std::cout << "u = " << u << "v = " << v << "\n";*/
      int j = 0;
      int max_edge_val = std::numeric_limits<int>::min();
      for(int i = 0; i < (n/2); i++){ 
        auto new_edge = rand() % 100, new_edge2 = rand() % 100; 
        tree.link(i, (2*i) + 1, new_edge);
        tree.link(i, (2*i) + 2, new_edge2);
        if(i == path[j]){
          if(i >= u && i < v){
            max_edge_val = std::max(max_edge_val, new_edge2);
          }
          j++;
        }
      }

      // Test return of max_edge_value.
      auto returned_query = tree.path_query(u, v);
      if(returned_query != max_edge_val){
        std::cout << seed << "\n";
        std::cout << "u = " << u << "v = " << v << "\n";
      }
      ASSERT_EQ(returned_query, max_edge_val);
    }
  }
}