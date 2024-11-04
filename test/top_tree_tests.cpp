#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/top_tree.h"

TEST(Top_tree_suite, constructor_test){
    TopTree<int> t(3, PATH, [] (int a, int b){return std::min(a,b);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    ASSERT_EQ(t.t.num_vertices, 3);
}

TEST(Top_tree_suite, basic_link_test){
    TopTree<int> t(10, PATH, [] (int a, int b){return std::min(a,b);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    for(int i = 1; i < 9; i++) t.link(i,i+1, i);
    for(int i = 1; i < 9 ; i++) {
        for(int j = i + 1; j < 9; j++){
            ASSERT_TRUE(t.connected(i,j));
        }
    }
}

TEST(Top_tree_suite, testInsertDeleteRakeLinkedList){
  std::vector<vertex_t> to_test({2, 10, 100, 1000});
  for(auto llist_size : to_test){
    TopTree<int> tree(llist_size);
    for(int i = 0; i < llist_size - 1; i++){
        tree.link(i , i + 1, i + 1); 
    }
    
    for(int i = 0; i < llist_size - 1; i++){
        for(int j = i + 1; j < llist_size - 1; j++){
            ASSERT_TRUE(tree.connected(i,j));
        }
    }
    /*for(int i = 0; i < llist_size - 1; i++){
      tree.cut(i , i + 1);
      ASSERT_FALSE(tree.connected(i, i+1));
    }*/
  }
}