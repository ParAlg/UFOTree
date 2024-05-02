#include <cstdlib>
#include <gtest/gtest.h>
#include "../include/rc_tree.h"
void create_tree1();
void print_tree(RCTree<int> *tree, int round);
//RCTree<int> *treeptr = nullptr;

void create_tree1() {
  RCTree<int> tree(7,3);
  RCCluster<int> edgeAB(9);
  edgeAB.boundary_vertexes.push_back(2);
  edgeAB.boundary_vertexes.push_back(0);
  

  RCCluster<int> edgeBC(1);
  edgeBC.boundary_vertexes.push_back(2);
  edgeBC.boundary_vertexes.push_back(1);

  RCCluster<int> edgeCE(2);
  edgeCE.boundary_vertexes.push_back(2);
  edgeCE.boundary_vertexes.push_back(4);
  
  RCCluster<int> edgeDE(3);
  edgeDE.boundary_vertexes.push_back(3);
  edgeDE.boundary_vertexes.push_back(4);

  RCCluster<int> edgeEF(6);
  edgeEF.boundary_vertexes.push_back(4);
  edgeEF.boundary_vertexes.push_back(5);

  RCCluster<int> edgeFG(5);
  edgeFG.boundary_vertexes.push_back(5);
  edgeFG.boundary_vertexes.push_back(6);

  tree.add_neighbor(0, &edgeAB, -1);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeBC, -1);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeCE, -1);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeDE, -1);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeEF, -1);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeFG, -1);
  //print_tree(tree, 0);  
}

void print_tree(RCTree<int> *tree, int round){
  std::cout << "Round: " << round << "\n";
  for(int i = 0; i < tree->n; i++){
    if(tree->contraction_tree[i].size() >= round + 1 && tree->contraction_tree[i][round] != nullptr){
      std::cout << "Vertex: " << i << " :[";
      for(int j = 0; j < tree->degree_bound; j++){
        if(tree->contraction_tree[i][round][j] != nullptr){
          //std::cout << "HERE!";
          auto cluster = tree->contraction_tree[i][round][j];
          std::cout << "[";
          for(int k = 0; k < cluster->boundary_vertexes.size(); k++){
            std::cout << cluster->boundary_vertexes[k] << ", "; 
          }
          std::cout << "]";
        }
      }
      std::cout << "] \n";
    }
  } 
}

TEST(RCTreeSuite, test_constructor) {
  RCTree<int> tree1(7, 3);
  ASSERT_EQ(tree1.degree_bound, 3);
  ASSERT_EQ(tree1.n, 7);
}

TEST(RCTreeSuite, test_helper_methods){
  RCTree<int> tree(7,3);
  RCCluster<int> edge1(9);
  RCCluster<int> unary1(13);
  RCCluster<int> unary2(14);
  edge1.boundary_vertexes.push_back(0);
  edge1.boundary_vertexes.push_back(1);

  unary1.boundary_vertexes.push_back(0);
  unary2.boundary_vertexes.push_back(0);

  //2 unary, 1 binary cluster. Degree = 1.
  tree.add_neighbor(0, &edge1, -1);
  tree.add_neighbor(0, &unary1, -1);
  tree.add_neighbor(0, &unary2, -1);

  //print_tree(tree, 0);
  ASSERT_EQ(tree.get_degree(0,0), 1);
  ASSERT_EQ(tree.get_degree(1,0), 1);
  //ASSERT_EQ(tree.neighbor_count(tree.adj[0][0]), 3);
  //ASSERT_EQ(tree.neighbor_count(tree.adj[0][1]), 1);
  ASSERT_EQ(tree.contracts(0, 0, 0), true);
  ASSERT_EQ(tree.contracts(0, 1, 0), false);
}

/*
TEST(RCTreeSuite, test_spread_affection){

}

TEST(RCTreeSuite, test_MIS){

}
*/
/*
TEST(RCTreeSuite, testRakeLinkedList){

  for(int llist_size = 1; llist_size < 100000; llist_size++){
    RCTree<int> tree(llist_size, 3);
    for(int i = 0; i < llist_size - 1; i++){
      tree.link(i, i+1, i+1);
    }
  }
}
*/

TEST(RCTreeSuite, testTHETREE){
  RCTree<int> tree(12, 3);
  tree.link(0,1,1);
  tree.link(1,2,2);
  tree.link(1,3,3);
  tree.link(3,4,4);
  tree.link(4,5,5);
  tree.link(4,7,6);
  tree.link(6,7,7);
  tree.link(7,8,8);

  tree.link(8,9,9);

  print_tree(&tree, 0);
  print_tree(&tree, 1);
  print_tree(&tree, 2);
  print_tree(&tree, 3);

  tree.link(8,10,10);
  tree.link(10,11,11);
}
