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