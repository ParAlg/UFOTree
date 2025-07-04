#include <gtest/gtest.h>
#include <unordered_set>
#include "../include/fast_ufo_tree.h"

using namespace dgbs;


template<typename v_t, typename e_t>
bool IUFOTree<v_t, e_t>::is_valid() {
}

template<typename v_t, typename e_t>
void IUFOTree<v_t, e_t>::print_tree() {
}

TEST(IUFOTreeSuite, cluster_constructor_test){
  IUFOCluster<int, int> c(1, 0, -1);
  ASSERT_TRUE(c.id == 1 && c.level == 0 && c.parent == -1);
}

TEST(IUFOTreeSuite, remove_ancestors_correctness_1){
  IUFOTree<int, int> tree(8);

  IUFOCluster<int, int> c1(0, 0, 1);
  IUFOCluster<int, int> c2(1, 1, 2);
  c2.tail = 0;
  
  IUFOCluster<int, int> c3(2, 2, 6);
  c3.tail = 1;

  IUFOCluster<int, int> c4(3, 0, 2);
  IUFOCluster<int, int> c5(4, 1, 6);
  IUFOCluster<int, int> c6(5, 0, 6);
  IUFOCluster<int, int> c7(6, 4, -1);
  c7.tail = 7;

  IUFOCluster<int, int> c8(7, 3, 6);
  
  // Cluster 3 linked list
  c2.prev = &c4;
  c4.next = &c2;

  // Cluster 7 linked list
  c8.prev = &c3;
  c3.next = &c8;
  c3.prev = &c5;
  c5.next = &c3;
  c5.prev = &c6;
  c6.next = &c5;
  
  tree.clusters[0] = &c1;
  tree.clusters[1] = &c2;
  tree.clusters[2] = &c3;
  tree.clusters[3] = &c4;
  tree.clusters[4] = &c5;
  tree.clusters[5] = &c6;
  tree.clusters[6] = &c7;
  tree.clusters[7] = &c8;

  tree.remove_ancestors(5,6);
  
  ASSERT_EQ(c1.parent, 1);
  ASSERT_EQ(c1.next, nullptr);
  ASSERT_EQ(c1.prev, nullptr);
  
  ASSERT_EQ(c2.parent, 2);
  ASSERT_EQ(c2.tail, 0);
  ASSERT_EQ(c2.next, nullptr);
  ASSERT_EQ(c2.prev, &c4);

  ASSERT_EQ(c3.parent, -1);
  ASSERT_EQ(c3.tail, 1);
  ASSERT_EQ(c3.next, nullptr);
  ASSERT_EQ(c3.prev, nullptr);

  ASSERT_EQ(c4.parent, 2);
  ASSERT_EQ(c4.tail, -1);
  ASSERT_EQ(c4.next, &c2);
  ASSERT_EQ(c4.prev, nullptr);

  ASSERT_EQ(c5.parent, -1);
  ASSERT_EQ(c5.level, 1);

  ASSERT_EQ(c6.parent, -1);
  ASSERT_EQ(c6.level, 0);

  ASSERT_EQ(c7.level, 0);
  ASSERT_EQ(c7.parent, -1);

  ASSERT_EQ(c8.level, 3);
  ASSERT_EQ(c8.parent, -1);
}
