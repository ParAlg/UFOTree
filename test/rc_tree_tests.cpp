#include <gtest/gtest.h>
#include "../include/rc_tree.h"


TEST(RCTreeSuite, test_constructor) {
  RCTree<int> tree1(7, 3);
  ASSERT_EQ(tree1.degree_bound, 3);
  ASSERT_EQ(tree1.n, 7);
}

TEST(RCTreeSuite, test_helper_methods){
  RCTree<int> tree(7, 3);
  RCCluster<int> edge1(9);
  RCCluster<int> unary1(13);
  RCCluster<int> unary2(14);
  edge1.boundary_vertexes.push_back(0);
  edge1.boundary_vertexes.push_back(1);

  unary1.boundary_vertexes.push_back(0);
  unary2.boundary_vertexes.push_back(0);
  
  //2 unary, 1 binary cluster. Degree = 1.
  tree.add_neighbor(tree.adj[0][0], &edge1);
  tree.add_neighbor(tree.adj[0][0], &unary1);
  tree.add_neighbor(tree.adj[0][0], &unary2);
  
  //1 binary cluster, nothing else. degree = 1
  tree.add_neighbor(tree.adj[0][1], &edge1);
  
  RCCluster<int>*** a = (RCCluster<int>***) calloc(7, sizeof(RCCluster<int>***));
  tree.adj.push_back(a);
  for(int i = 0; i < 7; i++){ tree.adj[1][i] = (RCCluster<int>**) calloc(3, sizeof(RCCluster<int>*));}
  ASSERT_EQ(tree.get_degree(0,0), 1);
  ASSERT_EQ(tree.get_degree(1,0), 1);
  ASSERT_EQ(tree.neighbor_count(tree.adj[0][0]), 3);
  ASSERT_EQ(tree.neighbor_count(tree.adj[0][1]), 1);
  ASSERT_EQ(tree.contracts(0, 0), true);
  ASSERT_EQ(tree.contracts(0, 1), false);
}

TEST(RCTreeSuite, test_spread){

}
