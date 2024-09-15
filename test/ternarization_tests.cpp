#include <gtest/gtest.h>
#include <ostream>
#include <stdexcept>
#include <random>
#include "../include/ternarized_tree.h"
#include "../include/rc_tree.h"
#include "../include/topology_tree.h"


template<typename DynamicTree, typename aug_t>
int TernarizedTree<DynamicTree, aug_t>::get_length_of_chain(vertex_t v){
  int length = 0;
  if(chain_map.find(v) == chain_map.end()) return length;
  auto curr = chain_map[v].second;
  while(curr != MAX_VERTEX_T){
    curr = chain_map[curr].second;
    length++;
  }
  return length;
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::vertex_on_chain(vertex_t start, vertex_t to_find){
  if(start > to_find) std::swap(start, to_find);
  if(edge_map.find(std::pair(start, to_find)) == edge_map.end()) return false;
  auto pair = edge_map[std::pair(start, to_find)];
  
  auto curr1 = chain_map.find(pair.first) != chain_map.end() && pair.first != start ? chain_map[pair.first].first : pair.first;
  auto curr2 = chain_map.find(pair.second) != chain_map.end() && pair.second != to_find ? chain_map[pair.second].first : pair.second;
  
  while(curr1 != start && chain_map[curr1].first != MAX_VERTEX_T){
    curr1 = chain_map[curr1].first;
  }
  while(curr2 != to_find && chain_map[curr2].first != MAX_VERTEX_T){
    curr2 = chain_map[curr2].first;
  }
  return curr1 == start && curr2 == to_find && tree.connected(start, to_find);
}

template<typename DynamicTree, typename aug_t>
bool TernarizedTree<DynamicTree, aug_t>::test_vertex_deleted(vertex_t u, vertex_t v){
  if(u > v) std::swap(u,v);
  if(edge_map.find(std::pair(u,v)) != edge_map.end()) return false;
  if(tree.connected(u,v)) return false;
  if(chain_map.find(u) != chain_map.end()){
    auto curr1 = u;
    while(chain_map[curr1].second != MAX_VERTEX_T){
      curr1 = chain_map[curr1].second;
      if(tree.connected(curr1, v)) return false;
    }
  }
  if(chain_map.find(v) != chain_map.end()){
    auto curr2 = v;
    while(chain_map[curr2].second != MAX_VERTEX_T){
      curr2 = chain_map[curr2].second;
      if(tree.connected(curr2, u)) return false;
    }
  }
  return true;
}
TEST(TernarizationSuite, constructor_test){
  
  TernarizedTree<RCTree<int>, int> tree(5);
  TernarizedTree<TopologyTree<int>, int> tree2(5);

  ASSERT_EQ(tree.n, 5);
  ASSERT_EQ(tree.tree.n, 10);
  ASSERT_EQ(tree2.n, 5);
}

TEST(TernarizationSuite, test_determine_link_v){ 
  // Case where vertex is unternarized.
  TernarizedTree<TopologyTree<int>, int> t(5);
  t.link(0, 1, 1);
  t.link(0, 2, 2);
  t.link(0, 3, 3);  
  ASSERT_TRUE(t.determine_link_v(0) == 6);
  
  // Case where vertex was already ternarized and now we have to add to head.
  TernarizedTree<TopologyTree<int>, int> t2(10);
  t2.link(0,1,1);
  t2.link(0,2,2);
  t2.link(0,3,3);
  t2.link(0,4,4); 
  ASSERT_EQ(t2.determine_link_v(0), 12);

  // Same cases as above but with RCTrees
  TernarizedTree<RCTree<int>, int> rt(5);
  rt.link(0, 1, 1);
  rt.link(0, 2, 2);
  rt.link(0, 3, 3);  
  ASSERT_TRUE(rt.determine_link_v(0) == 6);

  // Case where vertex was already ternarized and now we have to add to head.
  TernarizedTree<RCTree<int>, int> rt2(10);
  rt2.link(0,1,1);
  rt2.link(0,2,2);
  rt2.link(0,3,3);
  rt2.link(0,4,4);
  
  ASSERT_TRUE(rt2.determine_link_v(0) == 12);
}

// 3 different test functions for link, all test 2 possible cases to consider.
TEST(TernarizationSuite, test_link_cases_1){
  // Case 1, neither u nor v is degree 3, so normal link should take place.
  TernarizedTree<TopologyTree<int>, int> t(10);
  TernarizedTree<RCTree<int>, int> rt(10);

  t.link(0,1,1); rt.link(0,1,1);
  t.link(0,2,2); rt.link(0,2,2);
  t.link(4,3,4); rt.link(4,3,4);
  t.link(3,5,5); rt.link(3,5,5);
  t.link(0,3,3); rt.link(0,3,3);
  
  ASSERT_TRUE(t.edge_map[std::pair(0,1)].first == 0 && t.edge_map[std::pair(0,1)].second == 1);
  // Link between 1 degree 3 vertex and 1 non-degree 3 vertex.
  t.cut(0,3); rt.cut(0,3);
  t.link(0,6,6); rt.link(0,6,6); // Zero is now degree 3
  t.link(0,7,7); rt.link(0,7,7); // Zero degree 3 -> 4, 7 degree 1
  ASSERT_TRUE(t.vertex_on_chain(0,7));
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_TRUE(rt.vertex_on_chain(0,7));
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
}

TEST(TernarizationSuite, test_link_cases_2){
  // Link between 2 degree 3 vertices, (both should become connected via a new ternarized path.)
  TernarizedTree<TopologyTree<int>, int> t(20);
  TernarizedTree<RCTree<int>, int> rt(20);

  t.link(0,1,1); rt.link(0,1,1);
  t.link(0,2,2); rt.link(0,2,2);
  t.link(0,3,3); rt.link(0,3,3);

  t.link(4,5,5); rt.link(4,5,5);
  t.link(4,6,6); rt.link(4,6,6);
  t.link(4,7,7); rt.link(4,7,7);

  t.link(0,4,4); rt.link(0,4,4);

  ASSERT_TRUE(t.vertex_on_chain(0,4));
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_EQ(t.get_length_of_chain(4), 2);
  ASSERT_TRUE(rt.vertex_on_chain(0,4));
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
  ASSERT_EQ(rt.get_length_of_chain(4), 2);
 
  // Link between 2 previously ternarized vertices 
  t.cut(0,4); rt.cut(0,4);
  ASSERT_FALSE(t.vertex_on_chain(0,4));
  t.link(0,4,4); rt.link(0,4,4);
  ASSERT_TRUE(t.vertex_on_chain(0,4));
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_EQ(t.get_length_of_chain(4), 2);
  ASSERT_TRUE(rt.vertex_on_chain(0,4));
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
  ASSERT_EQ(rt.get_length_of_chain(4), 2);

  // Unternarized vertex becomes ternarized and is on the ternarized chain of another vertex
  t.link(1,8,8); rt.link(1,8,8);
  t.link(1,9,9); rt.link(1,9,9);
  t.link(1,10,10); rt.link(1,10,10);
  ASSERT_TRUE(t.vertex_on_chain(0,1));
  ASSERT_TRUE(rt.vertex_on_chain(0,1));
  ASSERT_EQ(t.get_length_of_chain(1), 2);
  ASSERT_EQ(rt.get_length_of_chain(1), 2);
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
}
TEST(TernarizationSuite, test_cut_cases_1){
  TernarizedTree<TopologyTree<int>, int> t(10);
  TernarizedTree<RCTree<int>, int> rt(10);
  // Cut between 2 unternarized vertices.
  t.link(0,1,1); rt.link(0,1,1);
  t.link(0,2,1); rt.link(0,2,1);
  t.link(0,3,1); rt.link(0,3,1);
  t.link(3,4,1); rt.link(3,4,1);
  t.link(3,5,1); rt.link(3,5,1);
  t.cut(0,3); rt.cut(0,3);
  ASSERT_TRUE(!t.connected(0,3) && !rt.connected(0,3));

  // Cut between vertex with a ternarized chain and one with no ternarized chain.
  t.link(0,6,6); rt.link(0,6,6);
  t.link(0,7,7); rt.link(0,7,7);
  t.link(0,3,3); rt.link(0,3,3);
  t.cut(0,3); rt.cut(0,3);
  ASSERT_FALSE(t.vertex_on_chain(0,3));
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_FALSE(t.connected(0,3));
  ASSERT_FALSE(rt.vertex_on_chain(0,3));
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
  ASSERT_EQ(rt.edge_map.find(std::pair(0,3)), rt.edge_map.end());
  ASSERT_EQ(t.edge_map.find(std::pair(0,3)), t.edge_map.end());
  
  // Cut between 2 vertices that are on each other's ternarized chain.
  t.link(6,8,8); rt.link(6,8,8);
  t.link(6,9,9); rt.link(6,9,9);
  t.link(6,3,3); rt.link(6,3,3);
  t.cut(0,6); rt.cut(0,6);
  ASSERT_FALSE(t.vertex_on_chain(0,6));
  ASSERT_FALSE(rt.vertex_on_chain(0,6));
  ASSERT_EQ(t.get_length_of_chain(0), 2);
  ASSERT_EQ(t.get_length_of_chain(6), 1);
  ASSERT_FALSE(t.connected(0,6));
  ASSERT_EQ(rt.get_length_of_chain(0), 2);
  ASSERT_EQ(rt.get_length_of_chain(6), 1);
  ASSERT_EQ(rt.edge_map.find(std::pair(0,6)), rt.edge_map.end());
  ASSERT_EQ(t.edge_map.find(std::pair(0,6)), t.edge_map.end());
}

TEST(TernarizationSuite, test_random_maximal_star_graphs){
  int num_trials = 1, max_n = 1024;
  int seeds[num_trials];
  srand(time(NULL));
  for(int i = 0; i < num_trials; i++){seeds[i] = rand();}
  for(int trial = 0; trial < num_trials; trial++){
    srand(seeds[trial]);
    std::mt19937 gen(seeds[trial]);
    auto t_size = rand() % max_n + 1;
    if(t_size < 10) continue;
    std::uniform_int_distribution<> distr((t_size/2), t_size);
    int max_degree = distr(gen); 
    if(max_degree == t_size) max_degree--;
    //std::cout << "Tsize, max_degree: " << t_size << " " << max_degree << "\n";
    std::vector<std::pair<vertex_t,vertex_t>> links;
    TernarizedTree<TopologyTree<int>, int> t(t_size);
    TernarizedTree<RCTree<int>, int> rt(t_size);
    for(int i = 0; i < t_size; i += (max_degree - 1)){
      if(i == max_degree - 1) i++;
      for(int j = i + 1; j < std::min(i + (max_degree - 1), t_size); j++){
        if(!t.connected(i,j)){
          try{
            t.link(i,j,j);
            rt.link(i,j,j);
            links.push_back(std::pair(i,j));
          } catch(std::invalid_argument()){
            std::cout << "Link failed on vertex " << i << " and " << j << "Seed, t_size, max_degree: " 
              << seeds[trial] << "," << t_size << "," << max_degree << "\n";
          }
        }
      }
    }
    t.link(0, max_degree, max_degree); rt.link(0,max_degree,max_degree);
    links.push_back(std::pair(0,max_degree));

    // Test to see if all vertices are connected via a chain.
    for(auto pair : links){
      ASSERT_TRUE(t.vertex_on_chain(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " not found in topology tree" <<", Seed: " << seeds[trial];
      ASSERT_TRUE(rt.vertex_on_chain(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " not found in RC tree" << ", Seed: " << seeds[trial];
    }

    // Delete all links.
    for(auto pair: links){
      try{
        t.cut(pair.first, pair.second); rt.cut(pair.first, pair.second);
      } catch(std::invalid_argument){
        std::cout << "Cut failed on " << pair.first << " " << pair.second << "Seed: " << seeds[trial];
      }
      ASSERT_TRUE(t.test_vertex_deleted(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " still connected in topology tree" << ", Seed: " << seeds[trial];
      ASSERT_TRUE(rt.test_vertex_deleted(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " still connected in RC tree" << ", Seed: " << seeds[trial];
    }

    // Make sure no vertices are still connected.
    for(int i = 0; i < t_size; i++){
      ASSERT_EQ(t.tree.get_degree(i), 0);
      ASSERT_EQ(t.get_length_of_chain(i), 0);
      ASSERT_EQ(rt.tree.get_degree(i), 0);
      ASSERT_EQ(rt.get_length_of_chain(i), 0);
    }
  }
}

TEST(TernarizationSuite, incremental_decremental_test_random_trees){
  // 1) Store all links in an array, make sure links result in vertices on chain
  // 2) Make sure length of chain is what is expected i.e. degree - 2
  int num_trials = 1, max_n = 1024;
  int seeds[num_trials];
  srand(time(NULL));
  for(int i = 0; i < num_trials; i++) seeds[i] = rand();
  for(int trial = 0; trial < num_trials; trial++){
    //std::cout << "Seed: " << seeds[trial] << "\n";
    srand(seeds[trial]);
    auto t_size = rand() % max_n;
    std::vector<std::pair<vertex_t,vertex_t>> links;
    TernarizedTree<TopologyTree<int>, int> t(t_size);
    TernarizedTree<RCTree<int>, int> rt(t_size);
    while(links.size() < t_size - 1){
      vertex_t u = rand() % t_size, v = rand() % t_size;
      if(u!= v && !t.connected(u,v)) {
        t.link(u,v,u); rt.link(u,v,u);
        links.push_back(std::pair(u,v));
      }
    }

    for(auto pair : links){
      ASSERT_TRUE(t.vertex_on_chain(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " not found in topology tree" <<", Seed: " << seeds[trial];
      ASSERT_TRUE(rt.vertex_on_chain(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " not found in RC tree" << ", Seed: " << seeds[trial];
    }

    while(links.size() > 0){ 
      int link_to_del = rand() % links.size();
      auto pair = links[link_to_del];
      t.cut(pair.first, pair.second); rt.cut(pair.first, pair.second);
      ASSERT_TRUE(t.test_vertex_deleted(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " still connected in topology tree" << ", Seed: " << seeds[trial];
      ASSERT_TRUE(rt.test_vertex_deleted(pair.first, pair.second)) 
        << "Vertices " << pair.first << " and " << pair.second << " still connected in RC tree" << ", Seed: " << seeds[trial];
      links.erase(links.begin() + link_to_del);
    }
    // Make sure no vertices are still connected.
    for(int i = 0; i < t_size; i++){
      ASSERT_EQ(t.tree.get_degree(i), 0);
      ASSERT_EQ(t.get_length_of_chain(i), 0);
      ASSERT_EQ(rt.tree.get_degree(i), 0);
      ASSERT_EQ(rt.get_length_of_chain(i), 0);
    }
  }
}

TEST (TernarizationSuite, kary_tree_correctness_test) {
  vertex_t n = 256;
  vertex_t k = 10;
  TernarizedTree<TopologyTree<int>, int> tt(n);
  TernarizedTree<RCTree<int>, int> rt(n);
  for (vertex_t i = 1; i < n; i++) {
    tt.link(i,(i-1)/k);
    for (vertex_t j = 1; j <= i; j++) ASSERT_TRUE(tt.connected(j,(j-1)/k));
    rt.link(i,(i-1)/k);
    for (vertex_t j = 1; j <= i; j++) ASSERT_TRUE(rt.connected(j,(j-1)/k));
  }
  for (vertex_t i = 1; i < n; i++) {
    tt.cut(i,(i-1)/k);
    for (vertex_t j = 1; j <= i; j++) ASSERT_TRUE(!tt.connected(j,(j-1)/k));
    for (vertex_t j = i+1; j < n; j++) ASSERT_TRUE(tt.connected(j,(j-1)/k));
    rt.cut(i,(i-1)/k);
    for (vertex_t j = 1; j <= i; j++) ASSERT_TRUE(!rt.connected(j,(j-1)/k));
    for (vertex_t j = i+1; j < n; j++) ASSERT_TRUE(rt.connected(j,(j-1)/k));
  }
}
