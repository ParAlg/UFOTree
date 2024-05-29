#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <unordered_set>
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
 * 1) all contracting neighbors are affected, degree 1, degree 2, degree 3
 * 2) Does not spread affection if all contracting are not affected
 * */
TEST(RCTreeSuite, test_spread_affection_contracting_affected){
  
  RCTree<int> tree(4, 3);

  RCCluster<int> e1(9);
  e1.boundary_vertexes = std::vector<int>({1,2});
  RCCluster<int> e2(2);
  e2.boundary_vertexes = std::vector<int>({1,3});
  
  //Degree 1
  tree.round_contracted[1] = 0;
  tree.round_contracted[2] = 1;
  tree.round_contracted[3] = 1;
  tree.add_neighbor(0, &e1, -1);
  tree.add_neighbor(0, &e2, -1);
  (tree.affected)->insert(1);
  tree.spread_affection(0);
   
  std::unordered_set<int> actual_affected({1,2,3});
  ASSERT_EQ(*(tree.affected), actual_affected);
  
  // Degree 2
  RCTree<int> tree2(4,3);
  tree2.add_neighbor(0, &e1, -1);
  tree2.add_neighbor(0, &e2, -1);
  (tree2.affected)->insert(2);
  (tree2.affected)->insert(3);
  tree2.round_contracted[2] = 0;
  tree2.round_contracted[3] = 0;
  tree2.round_contracted[1] = 1;
  tree2.spread_affection(0);
  std::unordered_set<int> actual_affected2({1,2,3});
  ASSERT_EQ(*(tree2.affected), actual_affected2);

  //Degree 3
  RCTree<int> tree3(5, 3);
  RCCluster<int> e3(2);
  e3.boundary_vertexes = std::vector<int>({1,0});
  tree3.add_neighbor(0, &e1, -1);
  tree3.add_neighbor(0, &e2, -1);
  tree3.add_neighbor(0, &e3, -1);
  (tree3.affected)->insert(2);
  (tree3.affected)->insert(3);
  (tree3.affected)->insert(0);
  tree3.round_contracted[0] = 0;
  tree3.round_contracted[2] = 0;
  tree3.round_contracted[3] = 0;
  tree3.round_contracted[1] = 1;
  tree3.spread_affection(0); 
  std::unordered_set<int> actual_affected3({1,2,3,0});
  ASSERT_EQ(*(tree3.affected), actual_affected3);
}

TEST(RCTreeSuite, spread_affection_uncontracting_not_affected){
  // Test to make sure no vertices that are not supposed to be affected 
  // are not getting affected.
  RCTree<int> tree(5,3);
 
  RCCluster<int> e1(9);
  e1.boundary_vertexes = std::vector<int>({1,2});
  RCCluster<int> e2(2);
  e2.boundary_vertexes = std::vector<int>({1,3});
  RCCluster<int> e3(3);
  e3.boundary_vertexes = std::vector<int>({0,1});

  //Degree 2
  tree.round_contracted[1] = 0;
  tree.round_contracted[2] = 1;
  tree.round_contracted[3] = 0;
  tree.add_neighbor(0, &e1, -1);
  (tree.affected)->insert(2);
  tree.spread_affection(0);
  
  std::unordered_set<int> actual_affected({2});
  ASSERT_EQ(*(tree.affected), actual_affected);
  
  //Degree 3
  tree.add_neighbor(0, &e2, -1);
  tree.add_neighbor(0, &e3, -1);
  tree.spread_affection(0);

  actual_affected.insert(3);
  tree.affected->insert(3);
  ASSERT_EQ(*(tree.affected), actual_affected);
}

/*
 * 1) No vertex already involved in a contraction gets selected
 * 2) Degree 3 vertices do not get selected
 * 3) Set is always maximal i.e. vertices that should've been selected are always selected
 * 4) Test the valid MIS function always returns a valid MIS
 */
TEST(RCTreeSuite, test_MIS){
  RCTree<int> linked_list(7, 3);
  
  RCCluster<int> e1(7);
  RCCluster<int> e2(6);
  RCCluster<int> e3(5);
  RCCluster<int> e4(4);
  RCCluster<int> e5(3);
  RCCluster<int> e6(2);
  
}

bool validTree(RCTree<int>* tree){
  auto t = *tree; 
  return true;
}
TEST(RCTreeSuite, testRakeLinkedList){
  std::unordered_set<int> to_test({1, 10, 100, 1000, 10000});
  
  for(auto llist_size : to_test){
    RCTree<int> tree(llist_size, 3);
    for(int i = 0; i < llist_size - 1; i++){
      try{
        tree.link(i , i + 1, i + 1);
      } catch(std::invalid_argument){
        auto height = tree.round_contracted[tree.root];
        for(int j = 0; j < height; j++){
          print_tree(&tree, j);
        }
        FAIL();
      }
    }
  }
}

TEST(RCTreeSuite, testTHETREE){
  RCTree<int> tree(12, 3);
  tree.link(0,1,1);
  tree.link(1,2,2);
  tree.link(1,3,3);
  tree.link(3,4,4);
  tree.link(4,5,5);
  tree.link(4,7,6); 
  tree.link(6,7,7);
  tree.cut(6,7);
  tree.link(7,8,8);
  tree.link(8,9,9);
  tree.link(8,10,10);
  tree.link(10,11,11);

}

TEST(RCTreeSuite, testInsertCompleteBinaryTree){
  std::unordered_set<int> to_test({1, 3, 31, 1023, 8191});
  for(auto n : to_test){
    RCTree<int> tree(n, 3);
    for(int i = 0; i < (n/2); i++){ 
      tree.link(i, (2*i) + 1, i);
      tree.link(i, (2*i) + 2, i);
    }
  }
}
