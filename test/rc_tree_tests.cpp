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

  tree.add_neighbor(0, &edgeAB);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeBC);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeCE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeDE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeEF);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeFG);
  //print_tree(tree, 0);  
}

void print_tree(RCTree<int> *treeptr, int round){
  RCTree<int> tree = *treeptr;
  std::cout << "Round: " << round << "\n";
  for(int i = 0; i < tree.n; i++){
    if(tree.adj[round][i] != nullptr){
      std::cout << "Vertex: " << i << " :[";
      for(int j = 0; j < tree.degree_bound; j++){
        if(tree.adj[round][i][j] != nullptr){
          //std::cout << "HERE!";
          auto cluster = *tree.adj[round][i][j];
          std::cout << "[";
          for(int k = 0; k < cluster.boundary_vertexes.size(); k++){
            std::cout << cluster.boundary_vertexes[k] << ", "; 
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
  tree.add_neighbor(0, &edge1);
  tree.add_neighbor(0, &unary1);
  tree.add_neighbor(0, &unary2);

  //print_tree(tree, 0);
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


TEST(RCTreeSuite, test_spread_affection){
  
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

  tree.add_neighbor(0, &edgeAB);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeBC);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeCE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeDE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeEF);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeFG);

  RCCluster<int>*** a = (RCCluster<int>***) calloc(tree.n, sizeof(RCCluster<int>***));
  tree.adj.push_back(a);
  for(int i = 0; i < 7; i++){ tree.adj[1][i] = (RCCluster<int>**) calloc(3, sizeof(RCCluster<int>*));}

  RCCluster<int> unaryA(9);
  unaryA.boundary_vertexes.push_back(2);

  RCCluster<int> unaryB(1);
  unaryB.boundary_vertexes.push_back(2);

  RCCluster<int> unaryD(3);
  unaryD.boundary_vertexes.push_back(4);

  RCCluster<int> binaryF(11);
  binaryF.boundary_vertexes.push_back(4);
  binaryF.boundary_vertexes.push_back(6);

  tree.add_neighbor(1, &unaryA);
  tree.add_neighbor(1, &unaryB);
  tree.add_neighbor(1, &edgeCE);
  tree.add_neighbor(1, &unaryD);
  tree.add_neighbor(1, &binaryF);
  
  print_tree(&tree, 0);
  print_tree(&tree, 1);

  tree.affected.insert(0);
  tree.affected.insert(1);
  tree.affected.insert(3);
  tree.affected.insert(5);

  tree.spread_affection(0);
  
  std::cout << "{";
  for(auto vertex : tree.affected){
    std::cout << vertex << " ";
  }
  std::cout << "}" << "\n";
}

TEST(RCTreeSuite, test_MIS){
 
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

  tree.add_neighbor(0, &edgeAB);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeBC);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeCE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeDE);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeEF);
  //print_tree(tree, 0);
  tree.add_neighbor(0, &edgeFG);
  
  RCCluster<int>*** a = (RCCluster<int>***) calloc(tree.n, sizeof(RCCluster<int>***));
  tree.adj.push_back(a);
  for(int i = 0; i < 7; i++){ tree.adj[1][i] = (RCCluster<int>**) calloc(3, sizeof(RCCluster<int>*));}

  RCCluster<int> unaryA(9);
  unaryA.boundary_vertexes.push_back(2);

  RCCluster<int> unaryB(1);
  unaryB.boundary_vertexes.push_back(2);

  RCCluster<int> unaryD(3);
  unaryD.boundary_vertexes.push_back(4);

  RCCluster<int> binaryF(11);
  binaryF.boundary_vertexes.push_back(4);
  binaryF.boundary_vertexes.push_back(6);

  tree.add_neighbor(1, &unaryA);
  tree.add_neighbor(1, &unaryB);
  tree.add_neighbor(1, &edgeCE);
  tree.add_neighbor(1, &unaryD);
  tree.add_neighbor(1, &binaryF);
  
  print_tree(&tree, 0);
  print_tree(&tree, 1);

  //tree.affected.insert(0);
  //tree.affected.insert(1);
  tree.affected.insert(2);
  tree.affected.insert(4);

  auto maximum_set = tree.MIS(0);
  for(auto v : maximum_set){
    std::cout << v << " ";
  }
}
