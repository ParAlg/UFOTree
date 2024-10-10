#include <gtest/gtest.h>
#include <limits>
#include <unordered_set>
#include "../include/topology_tree.h"


template<typename aug_t>
bool TopologyTree<aug_t>::is_valid() {
    std::unordered_set<vertex_t> cluster_set;
    std::unordered_set<vertex_t> next_cluster_set;
    for (vertex_t cluster_id = 0; cluster_id < clusters.size(); ++cluster_id) // Ensure that every pair of incident vertices are in the same component
        for (auto neighbor_id : clusters[cluster_id][0].neighbors) // This ensures all connectivity is correct by transitivity
            if (neighbor_id != NONE && get_root(cluster_id) != get_root(neighbor_id)) return false;
    for (int i = 0; i < clusters.size(); i++) cluster_set.insert(i);
    int level = 0;
    while (!cluster_set.empty()) {
        for (auto cluster_id : cluster_set) {
            auto cluster = &clusters[cluster_id][level];
            for (auto neighbor_id : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor_id != NONE && !clusters[neighbor_id][level].contains_neighbor(cluster_id)) return false;
            if (!contracts(cluster_id, level)) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    for (auto neighbor_id : cluster->neighbors)
                        if (neighbor_id != NONE && !contracts(neighbor_id, level)) return false;
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor_id : cluster->neighbors)
                        if (neighbor_id != NONE && !contracts(neighbor_id, level) && clusters[neighbor_id][level].get_degree() < 3) return false;
                } else if (cluster->get_degree() == 3) {
                    for (auto neighbor_id : cluster->neighbors)
                        if (neighbor_id != NONE && !contracts(neighbor_id, level) && clusters[neighbor_id][level].get_degree() < 2) return false;
                }
            }
            auto parent_id = get_parent(cluster_id, level);
            if (parent_id != NONE) next_cluster_set.insert(parent_id); // Get next level
        }
        cluster_set.swap(next_cluster_set);
        next_cluster_set.clear();
        level++;
    }
    return true;
}

template<typename aug_t>
int TopologyTree<aug_t>::get_height(vertex_t v) {
    int height = 0;
    vertex_t p = get_root(v);
    return clusters[p].size();
}

template<typename aug_t>
void TopologyTree<aug_t>::print_tree() {
    std::multimap<vertex_t, vertex_t> cluster_set;
    std::multimap<vertex_t, vertex_t> next_cluster_set;
    int level = 0;
    std::cout << "========================= LEAVES =========================" << std::endl;
    for (int i = 0; i < clusters.size(); i++) cluster_set.insert({get_parent(i, 0), i});
    for (auto entry : cluster_set) {
        auto leaf = entry.second;
        auto parent = entry.first;
        std::cout << "VERTEX " << leaf << " Parent " << parent << " Neighbors: ";
        for (auto neighbor_id : clusters[leaf][level].neighbors) if (neighbor_id != NONE) std::cout << neighbor_id << " ";
        std::cout << std::endl;
        bool in_map = false;
        for (auto entry : next_cluster_set) if (entry.second == parent) in_map = true;
        if (parent != NONE && !in_map) next_cluster_set.insert({get_parent(parent, level+1), parent});
    }
    cluster_set.swap(next_cluster_set);
    next_cluster_set.clear();
    level++;
    while (!cluster_set.empty()) {
        std::cout << "======================= NEXT LEVEL =======================" << std::endl;
        for (auto entry : cluster_set) {
            auto leaf = entry.second;
            auto parent = entry.first;
            std::cout << "VERTEX " << leaf << " Parent " << parent << " Neighbors: ";
            for (auto neighbor_id : clusters[leaf][level].neighbors) if (neighbor_id != NONE) std::cout << neighbor_id << " ";
            std::cout << std::endl;
            bool in_map = false;
            for (auto entry : next_cluster_set) if (entry.second == parent) in_map = true;
            if (parent != NONE && !in_map) next_cluster_set.insert({get_parent(parent, level+1), parent});
        }
        cluster_set.swap(next_cluster_set);
        next_cluster_set.clear();
        level++;
    }
}

TEST(TopologyTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int {return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int {return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(TopologyTreeSuite, incremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        TopologyTree<int> tree(n, qt, f, 0, 0);

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        std::vector<int> vertex_degrees(n,0);
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (vertex_degrees[u] >= 3) continue;
            if (vertex_degrees[v] >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << u << " and " << v << ".";
                links++;
                vertex_degrees[u]++;
                vertex_degrees[v]++;
            }
        }
    }
}

TEST(TopologyTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 128;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int {return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int {return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(TopologyTreeSuite, decremental_random_correctness_test) {
    int num_trials = 1;
    int seeds[num_trials];
    srand(time(NULL));
    for (int trial = 0; trial < num_trials; trial++) seeds[trial] = rand();
    for (int trial = 0; trial < num_trials; trial++) {
        vertex_t n = 256;
        QueryType qt = PATH;
        auto f = [](int x, int y)->int {return x + y;};
        TopologyTree<int> tree(n, qt, f, 0, 0);
        std::pair<vertex_t, vertex_t> edges[n-1];

        auto seed = seeds[trial];
        srand(seed);
        int links = 0;
        std::vector<int> vertex_degrees(n,0);
        while (links < n-1) {
            vertex_t u = rand() % n;
            vertex_t v = rand() % n;
            if (vertex_degrees[u] >= 3) continue;
            if (vertex_degrees[v] >= 3) continue;
            if (u != v && !tree.connected(u,v)) {
                tree.link(u,v);
                edges[links++] = {u,v};
                vertex_degrees[u]++;
                vertex_degrees[v]++;
            }
        }
        for (auto edge : edges) {
            auto u = edge.first;
            auto v = edge.second;
            tree.cut(u,v);
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
            ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << u << " and " << v << ".";
        }
    }
}

// TEST(TopologyTreeSuite, subtree_query_test) {
//     vertex_t n = 256;
//     QueryType qt = SUBTREE;
//     auto f = [](int x, int y)->int {return x + y;};
//     int id = 0;
//     int d = 1;
//     TopologyTree<int> tree(n, qt, f, id, d);

//     for (vertex_t i = 0; i < n-1; i++) {
//         tree.link(i,i+1);
//         for (vertex_t v = 0; v <= i+1; v++) {
//             ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i+2-v);
//             ASSERT_EQ(tree.subtree_query(v, (v <= i ? v+1 : MAX_VERTEX_T)), v+1);
//         }
//     }
//     for (vertex_t i = n-1; i > 0; i--) {
//         tree.cut(i-1,i);
//         for (vertex_t v = 0; v < i; v++) {
//             ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i-v);
//             ASSERT_EQ(tree.subtree_query(v, (v < i-1 ? v+1 : MAX_VERTEX_T)), v+1);
//         }
//     }
// }

// TEST(TopologyTreeSuite, path_query_test) {
//     vertex_t n = 64;
//     QueryType qt = PATH;
//     auto f = [](int x, int y)->int {return x + y;};
//     int id = 0;
//     int d = 0;
//     TopologyTree<int> tree(n, qt, f, id, d);

//     for (vertex_t i = 0; i < n-1; i++) {
//         tree.link(i,i+1,1);
//         for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
//             ASSERT_EQ(tree.path_query(u,v), v-u);
//     }
//     for (vertex_t i = 0; i < n-1; i++) {
//         tree.cut(i,i+1);
//         for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
//             ASSERT_EQ(tree.path_query(u,v), v-u);
//     }
// }

// TEST(TopologyTreeQuerySuite, LinkedListQueryTest){
//   std::vector<int> test_vals = {10, 100};
//   srand(time(NULL));
//   int seed = 1; 
//   srand(seed);
//   int num_trials = 100;
//   for(int n : test_vals){
//     for(int trial = 0; trial < num_trials; ++trial){
//       TopologyTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::min(x,y);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
//       //vertex_t u = rand() % (n-1), v = rand() % n; if(v == u){v++;}
//       vertex_t u = 4, v = 9;
//       if(v < u) std::swap(u,v);
//       /*std::cout << "Seed: " << seed << "\n";
//     std::cout << "u = " << u << "v = " << v << "\n";*/
//       int min_edge_val = std::numeric_limits<int>::max();
//       for(int i = 0; i < n-1; i++){
//         int new_edge = rand() % 100;
//         tree.link(i, i+1, new_edge);
//         if(i >= u && i < v) min_edge_val = std::min(min_edge_val, new_edge);
//       }

//       // Test return of min_edge_value.
//       auto returned_query = tree.path_query(u, v);
//       if(returned_query != min_edge_val){
//         std::cout << seed << "\n";
//         std::cout << "u = " << u << "v = " << v << "\n";
//       }
//       ASSERT_EQ(returned_query, min_edge_val);
//     }
//   }
// }

// TEST(TopologyTreeQuerySuite, BinaryTreeQueryTest){
//   std::vector<int> test_vals = {7, 31, 127};
//   srand(time(NULL));
//   int seed = 1; 
//   srand(seed);
//   int num_trials = 100;
//   for(int n : test_vals){
//     for(int trial = 0; trial < num_trials; ++trial){
//       TopologyTree<int> tree(n, QueryType::PATH, [] (int x, int y){return std::min(x,y);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
//       std::vector<int> path; 
//       for(int i = 0; i < n - 1; i = (2*i) + 2){
//         path.push_back(i);
//       }

//       vertex_t upper = rand() % (path.size() - 1), lower = rand() % path.size();if(lower == upper) lower++;
//       if(lower < upper) std::swap(lower, upper);

//       auto u = path[upper], v = path[lower];
//       //vertex_t u = 4, v = 7;
//       //if(v < u) std::swap(u,v);
//       /*std::cout << "Seed: " << seed << "\n";
//       std::cout << "u = " << u << "v = " << v << "\n";*/ 
//       int j = 0;
//       int min_edge_val = std::numeric_limits<int>::max();
//       for(int i = 0; i < (n/2); i++){ 
//         auto new_edge = rand() % 100, new_edge2 = rand() % 100; 
//         tree.link(i, (2*i) + 1, new_edge);
//         tree.link(i, (2*i) + 2, new_edge2);
//         if(i == path[j]){
//           if(i >= u && i < v){
//             min_edge_val = std::min(min_edge_val, new_edge2);
//           }
//           j++;
//         }
//       }

//       // Test return of min_edge_value.
//       auto returned_query = tree.path_query(u, v);
//       if(returned_query != min_edge_val){
//         std::cout << seed << "\n";
//         std::cout << "u = " << u << "v = " << v << "\n";
//       }
//       ASSERT_EQ(returned_query, min_edge_val);
//     }
//   }
// }
