#include <gtest/gtest.h>
#include "../include/ternarized_tree.h"
#include "../include/rc_tree.h"
#include "../include/topology_tree.h"

template<typename DynamicTree, typename TreeCluster, typename aug_t>
bool TernarizedTree<DynamicTree, TreeCluster, aug_t>::is_valid_ternarized_tree(){
  // To write function:
  // 1) Call is_valid_tree methods from Topology Tree and RCTree classes.
  // 2) Verify high degree vertices are ternarized.
  // 3) Verify length of ternarized chain.
  // 4) Verify chain contains vertex we expect it to.
}

template<typename DynamicTree, typename TreeCluster, typename aug_t>
int TernarizedTree<DynamicTree, TreeCluster, aug_t>::get_length_of_chain(vertex_t v){

}

template<typename DynamicTree, typename TreeCluster, typename aug_t>
bool TernarizedTree<DynamicTree, TreeCluster, aug_t>::vertex_on_chain(vertex_t start, vertex_t to_find){

}

TEST(TernarizationSuite, constructor_test){
  
  TernarizedTree<RCTree<int>, RCCluster<int>, int> tree(5);
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> tree2(5);

  ASSERT_EQ(tree.n, 5);
  ASSERT_EQ(tree.tree.n, 10);
  ASSERT_EQ(tree2.n, 5);
}

TEST(TernarizationSuite, test_get_id){
  //Topology Tree and RCCluster version of testing.
  TopologyCluster<int> cluster1(20, 3);
  RCCluster<int> cluster2(20);
  cluster2.boundary_vertexes.push_back(3); 
  cluster2.boundary_vertexes.push_back(4); 

  TopologyTree<int> t(2); 
  RCTree<int> t2(3);

  ASSERT_EQ(t.get_vertex_id(&cluster1, 3), 3);
  ASSERT_EQ(t2.get_vertex_id(&cluster2, 4), 3);
}

TEST(TernarizationSuite, test_get_edge_val){
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> t(5);
  t.link(0,1, 1); 
  t.link(0,2, 2);
  t.link(0,3, 3);
  ASSERT_EQ(t.get_edge_val(0), 1); 

  TernarizedTree<RCTree<int>, RCCluster<int>, int> rt(5);
  rt.link(0,1,1);
  rt.link(0,2,2);
  rt.link(0,3,3);
  ASSERT_EQ(rt.get_edge_val(0), 1);
}

TEST(TernarizationSuite, test_determine_link_v){ 
  // Case where vertex is unternarized.
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> t(5);
  t.link(0, 1, 1);
  t.link(0, 2, 2);
  t.link(0, 3, 3);  
  ASSERT_TRUE(t.determine_link_v(0) == 6);
  
  // Case where vertex was already ternarized and now we have to add to head.
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> t2(10);
  t2.link(0,1,1);
  t2.link(0,2,2);
  t2.link(0,3,3);
  t2.link(0,4,4); 
  ASSERT_TRUE(t2.determine_link_v(0) == 12);

  // Same cases as above but with RCTrees
  TernarizedTree<RCTree<int>, RCCluster<int>, int> rt(5);
  rt.link(0, 1, 1);
  rt.link(0, 2, 2);
  rt.link(0, 3, 3);  
  ASSERT_TRUE(rt.determine_link_v(0) == 6);

  // Case where vertex was already ternarized and now we have to add to head.
  TernarizedTree<RCTree<int>, RCCluster<int>, int> rt2(10);
  rt2.link(0,1,1);
  rt2.link(0,2,2);
  rt2.link(0,3,3);
  rt2.link(0,4,4);
  
  ASSERT_TRUE(rt2.determine_link_v(0) == 12);
}

// 3 different test functions for link, all test 2 possible cases to consider.
TEST(TernarizationSuite, test_link_cases_1){
  // Case 1, neither u nor v is degree 3, so normal link should take place.
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> t(10);
  TernarizedTree<RCTree<int>, RCCluster<int>, int> rt(10);

  t.link(0,1,1); rt.link(0,1,1);
  t.link(0,2,2); rt.link(0,2,2);
  t.link(4,3,4); rt.link(4,3,4);
  t.link(3,5,5); rt.link(3,5,5);
  t.link(0,3,3); rt.link(0,3,3);

  auto zero_neighbors = t.tree.get_neighbors(0);
  auto three_neighbors = t.tree.get_neighbors(3);
  auto zero_neighbors_rctree = rt.tree.get_neighbors(0);
  auto three_neighbors_rctree = rt.tree.get_neighbors(3);

  for(int i = 0; i < 3; i++){
    ASSERT_TRUE(t.get_id(zero_neighbors[i], 0) < 10);
    ASSERT_TRUE(t.get_id(three_neighbors[i], 3) < 10);
    ASSERT_TRUE(rt.get_id(zero_neighbors_rctree[i], 0) < 10);
    ASSERT_TRUE(rt.get_id(three_neighbors_rctree[i], 3) < 10);
  }  

  // Link between 1 degree 3 vertex and 1 non-degree 3 vertex.
  t.cut(0,3); rt.cut(0,3);
  t.link(0,6,6); rt.link(0,6,6); // Zero is now degree 3
  t.link(0,7,7); rt.link(0,7,7); // Zero degree 3 -> 4, 7 degree 1
  bool ternarized_v_found_top = false, ternarized_v_found_rc = false;
  for(int i = 0; i < 3; i++){
    if(t.get_id(zero_neighbors[i], 0) >= 10) ternarized_v_found_top = true;
    if(rt.get_id(zero_neighbors_rctree[i], 0) >= 10) ternarized_v_found_rc = true;
  }
  ASSERT_TRUE(ternarized_v_found_top && ternarized_v_found_rc);
}

TEST(TernarizationSuite, test_link_cases_2){
  // Link between 2 degree 3 vertices, (both should become connected via a new ternarized path.)
  TernarizedTree<TopologyTree<int>, TopologyCluster<int>, int> t(20);
  TernarizedTree<RCTree<int>, RCCluster<int>, int> rt(20);

  t.link(0,1,1); rt.link(0,1,1);
  t.link(0,2,2); rt.link(0,2,2);
  t.link(0,3,3); rt.link(0,3,3);

  t.link(4,5,5); rt.link(4,5,5);
  t.link(4,6,6); rt.link(4,6,6);
  t.link(4,7,7); rt.link(4,7,7);

  t.link(0,4,4); rt.link(0,4,4);

  ASSERT_TRUE(t.connected(0,4) && rt.connected(0,4));
  ASSERT_TRUE(t.edge_map[std::pair(0,4)].first >= 20 && t.edge_map[std::pair(0,4)].second >= 20);
  ASSERT_TRUE(rt.edge_map[std::pair(0,4)].first >= 20 && rt.edge_map[std::pair(0,4)].second >= 20);

  // Link between 2 previously ternarized vertices
  
}
TEST(TernarizationSuite, test_delete_ternarized_vertex){

}

TEST(TernarizationSuite, test_maximal_star_graphs){

}
