#include <cstdlib>
#include <gtest/gtest.h>
#include <stdexcept>
#include <unordered_set>
#include "../include/rc_tree.h"
void create_tree1();
void print_tree(RCTree<int> *tree, int round);
//RCTree<int> *treeptr = nullptr;

void create_tree1() {
  RCTree<int> tree(7);
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

template <typename aug_t>
void RCTree<aug_t>::is_valid_induced_tree(int round){
  /* This method verifies whether the induced tree constructed at a particular round 
   * was computed correctly. It checks for the following conditions: 
    Same edge has not been added twice for each vertex.
    Cluster has adjacency list its in as a boundary vertex, validating clusters.
    Comes from something at the previous level - i.e. if it is a contracted cluster, it was either present
                                                 at the level before or is from a vertex that contracted in
                                                 the previous round.
    Edges are bidirectional.*/

  for(int v = 0; v < n; v++){
    if(contraction_tree[v].size() >= round + 1){
      auto adjacencies = contraction_tree[v][round];
      for(int j = 0; j < degree_bound; j++){
        if(adjacencies[j] == nullptr) continue;
        auto deref_cluster = *adjacencies[j];
        auto boundaries = deref_cluster.boundary_vertexes;
        if(boundaries.size() == 1){
          if(boundaries[0] != v){
            throw std::invalid_argument("Unary cluster does not have correct boundary_v");
          }
        } else {
          if(boundaries[0] != v && boundaries[1] != v) 
            throw std::invalid_argument("Binary cluster does not have correct verts as boundary_vertexes");

          auto neighbor = GET_NEIGHBOR(v, deref_cluster);
          if(contraction_tree[neighbor].size() <= round) 
            throw std::invalid_argument("Neighbor did not live long enough.");

          bool contains_cluster = false;
          for(int i = 0; i < degree_bound; i++){
            if(contraction_tree[neighbor][round][i] == adjacencies[j]){
              contains_cluster = true; 
              break;
            }
          }
          if(!contains_cluster){
            throw std::invalid_argument("Unidirectional edge found, in valid_induced_tree");
          }
        }

        for(int i = 0; i < degree_bound; i++){
          if(i == j) continue;
          if(adjacencies[i] == adjacencies[j]){
            throw std::invalid_argument("Same cluster has been added twice");
          }
        }
        
        if(round > 0){
          bool in_prev_round = false;
          for(int i = 0; i < degree_bound; i++){
            if(contraction_tree[v][round - 1][i] == adjacencies[j]){
              in_prev_round = true;
              break;
            }

            if(contraction_tree[v][round - 1][i] != nullptr){
              auto neighbor_boundaries = contraction_tree[v][round - 1][i]->boundary_vertexes;
              if(std::find(neighbor_boundaries.begin(), neighbor_boundaries.end(), deref_cluster.rep_vertex)
                != neighbor_boundaries.end()){
                in_prev_round = true;
                break;
              }
            }
          }
          if(!in_prev_round){
            throw std::invalid_argument("Cluster comes from something that should not exist.");
          }
        }
      }
    }
  }
}

template<typename aug_t>
void RCTree<aug_t>::is_valid_MIS(int round){
  //Given a computed MIS, check to see if the MIS is maximal and no
  //vertices added are neighbors of each other.

  //Checks to see if the set is maximal by going through all vertexes in the
  //affected set and checking to see if atleast one neighbor of every vertex is in
  //the maximal independent set.
  for(vertex_t vertex : affected){
    // Each vertex being considered and its neighbors must survive to the present round
    // if an MIS is being computed, if this is not true then there is an error.
    if(contraction_tree[vertex][round] == nullptr){
      throw std::invalid_argument("Vertex " + std::to_string(vertex) + 
                                  " did not live at round " + std::to_string(round));
    }
    bool neighbor_in_max = false;
    for(int i = 0; i < degree_bound; i++){
      if(is_edge(contraction_tree[vertex][round][i])){
        auto dereferenced_cluster = *contraction_tree[vertex][round][i];
        int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
        if(maximal_set.count(vertex) == 0 && 
          ((contracts(neighbor, round, 0)  && affected.count(neighbor) == 0)|| 
          maximal_set.count(neighbor) == 1)){
          neighbor_in_max = true; // If either the vertex is next to an contracting unaffected neighbor or 
          // neighbor included in the MIS then increment count as this vertex is valid.
        }

        // If 2 adjacent neighbors have been added then also throw an exception.
        // as the set is not independent.
        if(maximal_set.count(vertex) == 1 && 
            maximal_set.count(neighbor) != 0){
          throw std::invalid_argument("2 adjacent vertices have been added to the MIS");
        }
      }
    }

    //If no neighbors are in the maximal set then this vertex should have been
    //in the maximal set. Thus the MIS is invalid
    if(maximal_set.count(vertex) == 0 && get_degree(vertex, round) < 3 &&
      !neighbor_in_max){
      throw std::invalid_argument("The set is not maximal");
    }
  }
}

template<typename aug_t>
void print_cluster(RCCluster<aug_t>* cluster){
  // Prints a cluster.
  std::cout << "Cluster: " << cluster->rep_vertex << "\n";
  std::cout << "Parent:" << cluster->parent << "\n";
  std::cout << "Aug_val : " << cluster->aug_val << "\n";
  std::cout << "boundary_vertexes: ";
  
  for(int i = 0; i < cluster->boundary_vertexes.size(); i++){
    std::cout << cluster->boundary_vertexes[i] << " ";
  } 
}
TEST(RCTreeSuite, test_constructor) {
  RCTree<int> tree1(7);
  ASSERT_EQ(tree1.degree_bound, 3);
  ASSERT_EQ(tree1.n, 7);
}

TEST(RCTreeSuite, test_helper_methods){
  RCTree<int> tree(7);
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
 * 
*/

TEST(RCTreeSuite, test_spread_affection_contracting_affected){

  RCTree<int> tree(4);

  RCCluster<int> e1(9);
  e1.boundary_vertexes = std::vector<vertex_t>({1,2});
  RCCluster<int> e2(2);
  e2.boundary_vertexes = std::vector<vertex_t>({1,3});

  //Degree 1
  tree.round_contracted[1] = 0;
  tree.round_contracted[2] = 1;
  tree.round_contracted[3] = 1;
  tree.add_neighbor(0, &e1, -1);
  tree.add_neighbor(0, &e2, -1);
  (tree.affected).insert(1);
  tree.spread_affection(0);

  std::unordered_set<vertex_t> actual_affected({1,2,3});
  ASSERT_EQ((tree.affected), actual_affected);

  // Degree 2
  RCTree<int> tree2(4);
  tree2.add_neighbor(0, &e1, -1);
  tree2.add_neighbor(0, &e2, -1);
  (tree2.affected).insert(2);
  (tree2.affected).insert(3);
  tree2.round_contracted[2] = 0;
  tree2.round_contracted[3] = 0;
  tree2.round_contracted[1] = 1;
  tree2.spread_affection(0);
  std::unordered_set<vertex_t> actual_affected2({1,2,3});
  ASSERT_EQ((tree2.affected), actual_affected2);

  //Degree 3
  RCTree<int> tree3(5);
  RCCluster<int> e3(2);
  e3.boundary_vertexes = std::vector<vertex_t>({1,0});
  tree3.add_neighbor(0, &e1, -1);
  tree3.add_neighbor(0, &e2, -1);
  tree3.add_neighbor(0, &e3, -1);
  (tree3.affected).insert(2);
  (tree3.affected).insert(3);
  (tree3.affected).insert(0);
  tree3.round_contracted[0] = 0;
  tree3.round_contracted[2] = 0;
  tree3.round_contracted[3] = 0;
  tree3.round_contracted[1] = 1;
  tree3.spread_affection(0); 
  std::unordered_set<vertex_t> actual_affected3({1,2,3,0});
  ASSERT_EQ((tree3.affected), actual_affected3);
}

TEST(RCTreeSuite, spread_affection_uncontracting_not_affected){
  // Test to make sure no vertices that are not supposed to be affected 
  // are not getting affected.
  RCTree<int> tree(5);

  RCCluster<int> e1(9);
  e1.boundary_vertexes = std::vector<vertex_t>({1,2});
  RCCluster<int> e2(2);
  e2.boundary_vertexes = std::vector<vertex_t>({1,3});
  RCCluster<int> e3(3);
  e3.boundary_vertexes = std::vector<vertex_t>({0,1});

  //Degree 2
  tree.round_contracted[1] = 0;
  tree.round_contracted[2] = 1;
  tree.round_contracted[3] = 0;
  tree.add_neighbor(0, &e1, -1);
  (tree.affected).insert(2);
  tree.spread_affection(0);

  std::unordered_set<vertex_t> actual_affected({2});
  ASSERT_EQ((tree.affected), actual_affected);

  //Degree 3
  tree.add_neighbor(0, &e2, -1);
  tree.add_neighbor(0, &e3, -1);
  tree.spread_affection(0);

  actual_affected.insert(3);
  tree.affected.insert(3);
  ASSERT_EQ((tree.affected), actual_affected);
}

/*
 * 1) No vertex already involved in a contraction gets selected
 * 2) Degree 3 vertices do not get selected
 * 3) Set is always maximal i.e. vertices that should've been selected are always selected
 * 4) Test the valid MIS function always returns a valid MIS
 */
TEST(RCTreeSuite, test_MIS){
  RCTree<int> linked_list(7);

  RCCluster<int> e1(7);
  RCCluster<int> e2(6);
  RCCluster<int> e3(5);
  RCCluster<int> e4(4);
  RCCluster<int> e5(3);
  RCCluster<int> e6(2);

}

TEST(RCTreeSuite, testBadCaseStarGraph){
  RCTree<int> tree(20);
  tree.link(0,1,1);
  tree.link(0,6,6);
  tree.link(0,11,11);
  for(int i = 1; i <= 4; i++){
    tree.link(1 + (i - 1), 1 + i, i);
    tree.link(6 + (i - 1), 6 + i, 6 + i);
    tree.link(11 + (i - 1), 11 + (i), 11 + i);
  }

  tree.cut(1,2);
}

TEST(RCTreeSuite, testInsertDeleteRakeLinkedList){
  std::unordered_set<vertex_t> to_test({2, 10, 100, 1000});
  for(auto llist_size : to_test){
    RCTree<int> tree(llist_size);
    for(int i = 0; i < llist_size - 1; i++){
      try{
        tree.link(i , i + 1, i + 1);
      } catch(std::invalid_argument){
        auto height = tree.round_contracted[tree.roots[0]];
        for(int j = 0; j < height; j++){
          print_tree(&tree, j);
        }
        FAIL();
      }
    }

    ASSERT_EQ(tree.representative_clusters[tree.roots[0]]->aug_val, (llist_size * (llist_size - 1))/2);
    for(int i = 0; i < llist_size - 1; i++){
      tree.cut(i , i + 1);
    }
  }
}

TEST(RCTreeSuite, testTHETREE){
  RCTree<int> tree(12);
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

TEST(RCTreeSuite, testInsertDeleteCompleteBinaryTree){
  std::unordered_set<vertex_t> to_test({1, 3, 31, 255});
  for(auto n : to_test){
    RCTree<int> tree(n);
    for(int i = 0; i < (n/2); i++){ 
      tree.link(i, (2*i) + 1, i);
      tree.link(i, (2*i) + 2, i);
    }

    for(int i = 0; i < (n/2); i++){
      tree.cut(i, (2*i) + 1);
      tree.cut(i, (2*i) + 2);
    }
  }
}

TEST(RCTreeSuite, randomDecrementalTestBinaryTree){
  int n = 1023; 
  RCTree<int> tree(n);
  for(int i = 0; i < (n/2); i++){ 
    tree.link(i, (2*i) + 1, i);
    tree.link(i, (2*i) + 2, i);
  }  

  for(int i = 0; i < 1023; i++){
    int u = std::rand() % 1023;
    if(tree.edge_exists(u, (2*u) + 1)) tree.cut(u, (2*u)+1);

    if(tree.edge_exists(u, (2*u) + 2)) tree.cut(u, (2*u) + 2);
  }
}

TEST(RCTreeSuite, randomIncrementalTestBinaryTree){
  int n = 1024;
  RCTree<int> tree(n);
  for(int i = 0; i < 10000; i++){
    int u = std::rand() % n;
    int v = std::rand() % n;

    if(!tree.connected(u, v) && tree.get_degree(u, 0) < 3 && tree.get_degree(v, 0) < 3) tree.link(u, v, u);
  }

  for(int i = 0; i < 10000; i++){
    int u = std::rand() % n;
    for(int j = 0; j < 3; j++){
      if(tree.is_edge(tree.contraction_tree[u][0][j])){
        auto dereferenced_cluster = *tree.contraction_tree[u][0][j];
        auto neighbor = GET_NEIGHBOR(u, dereferenced_cluster);
        tree.cut(u, neighbor);
      }
    }
  }
}

TEST(RCTreeSuite, randomIncrementalTests){
  int num_trials = 1;
  int seeds[num_trials];
  srand(time(NULL));
  for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
  for (int trial = 0; trial < num_trials; trial++) {
    vertex_t n = 256;
    RCTree<int> tree(n);

    auto seed = seeds[trial];
    srand(seed);
    int links = 0;
    while (links < n-1) {
      vertex_t u = rand() % n;
      vertex_t v = rand() % n;
      if (u != v && !tree.connected(u,v) && tree.get_degree(u, 0) < 3 && tree.get_degree(v, 0) < 3) {
        try {
          tree.link(u,v, u);
          links++;
        } catch(std::invalid_argument){
          std::cout << "Tree invalid after linking " << u << " and " << v << ".\n";
          print_tree(&tree, 0);
          break;
        }
      }
    }
  }
}

TEST(RCTreeSuite, decremental_random_correctness_test) {
  int num_trials = 1;
  int seeds[num_trials];
  srand(time(NULL));
  for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
  for (int trial = 0; trial < num_trials; trial++) {
    vertex_t n = 256;
    RCTree<int> tree(n);
    std::pair<vertex_t, vertex_t> edges[n-1];

    auto seed = seeds[trial];
    srand(seed);
    int links = 0;
    while (links < n-1) {
      vertex_t u = rand() % n;
      vertex_t v = rand() % n;
      if (u != v && !tree.connected(u,v) && tree.get_degree(u, 0) < 3 && tree.get_degree(v, 0) < 3) {
        tree.link(u,v,u);
        edges[links++] = {u,v};
      }
    }
    for (auto edge : edges) {
      auto u = edge.first;
      auto v = edge.second;
      try{
        tree.cut(u,v);
        ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
      } catch(std::invalid_argument){
        std::cout << "Tree invalid after cutting " << u << " and " << v << ".\n";
        print_tree(&tree, 0);
        break;
      }
    }
  }
}
