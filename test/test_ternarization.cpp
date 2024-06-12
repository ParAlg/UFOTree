#include <gtest/gtest.h>
#include "../include/ternarized_tree.h"
#include "../include/rc_tree.h"

template<typename DynamicTree, typename TreeCluster, typename aug_t>
bool TernarizedTree<DynamicTree, TreeCluster, aug_t>::is_valid_ternarized_tree(){
  
}
TEST(TernarizationSuite, constructor_test){
  
  TernarizedTree<RCTree<int>, RCCluster<int>, int> tree(5);
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> tree2(5);

  ASSERT_EQ(tree.n, 5);
  ASSERT_EQ(tree.tree.n, 10);
  ASSERT_EQ(tree2.n, 5);
}

TEST(TernarizationSuite, test_get_id){

}

TEST(TernarizationSuite, test_get_edge_val){

}

TEST(TernarizationSuite, test_determine_link_v){

}

TEST(TernarizationSuite, test_ternarize_vertex){

}

TEST(TernarizationSuite, test_delete_ternarized_vertex){

}

TEST(TernarizationSuite, test_maximal_star_graphs){

}
