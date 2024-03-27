#include <gtest/gtest.h>
#include "../include/topology_tree.h"


template<typename aug_t>
bool TopologyTree<aug_t>::is_valid() {
    std::unordered_set<TopologyCluster<aug_t>*> clusters;
    std::unordered_set<TopologyCluster<aug_t>*> next_clusters;
    for (int i = 0; i < this->leaves.size(); i++) clusters.insert(&this->leaves[i]);
    while (!clusters.empty()) {
        for (auto cluster : clusters) {
            for (auto neighbor : cluster->neighbors) // Ensure all neighbors also point back
                if (neighbor && !neighbor->contains_neighbor(cluster)) return false;
            if (!cluster->contracts()) { // Ensure maximality of contraction
                if (cluster->get_degree() == 1) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts()) return false;
                } else if (cluster->get_degree() == 2) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 3) return false;
                } else if (cluster->get_degree() == 3) {
                    for (auto neighbor : cluster->neighbors)
                        if (neighbor && !neighbor->contracts() && neighbor->get_degree() < 2) return false;
                }
            }
            if (cluster->parent) next_clusters.insert(cluster->parent); // Get next level
        }
        // Maximality should ensure this, but we leave the test as a sanity check
        if (6*next_clusters.size() > 5*clusters.size()) return false;
        clusters.swap(next_clusters);
        next_clusters.clear();
    }
    return true;
}

TEST(TopologyTreeSuite, incremental_binarytree_speed_test) {
    vertex_t n = 1000;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
}

TEST(TopologyTreeSuite, incremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, incremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        for (vertex_t u = 0; u < 2*i+2; u++) for (vertex_t v = u+1; v <= 2*i+2; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after linking " << i << " with children.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t u = 0; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
        ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all links.";
}

TEST(TopologyTreeSuite, decremental_linkedlist_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
        for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " and " << i+1 << ".";
    }
}

TEST(TopologyTreeSuite, decremental_binarytree_correctness_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f, 0, 0);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.cut(i,2*i+1);
        tree.cut(i,2*i+2);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
        ASSERT_FALSE(tree.connected(2*i+1,2*i+2)) << "Vertex " << 2*i+1 << " and " << 2*i+2 << " connected.";
        if (4*i+6 < n) {
            ASSERT_TRUE(tree.connected(4*i+3,4*i+4)) << "Vertex " << 4*i+3 << " and " << 4*i+4 << " not connected.";
            ASSERT_TRUE(tree.connected(4*i+5,4*i+6)) << "Vertex " << 4*i+3 << " and " << 4*i+4 << " not connected.";
        }
        ASSERT_TRUE(tree.is_valid()) << "Tree invalid after cutting " << i << " from children.";
    }
    if (n%2 == 0) tree.cut((n-1)/2,n-1);
    for (vertex_t u = 0; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
        ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
    ASSERT_TRUE(tree.is_valid()) << "Tree invalid after all cuts.";
}

TEST(TopologyTreeSuite, subtree_query_test) {
    vertex_t n = 256;
    QueryType qt = SUBTREE;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 1;
    TopologyTree<int> tree(n, qt, f, id, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t v = 0; v <= i+1; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i+2-v);
            ASSERT_EQ(tree.subtree_query(v, (v <= i ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
    for (vertex_t i = n-1; i > 0; i--) {
        tree.cut(i-1,i);
        for (vertex_t v = 0; v < i; v++) {
            ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i-v);
            ASSERT_EQ(tree.subtree_query(v, (v < i-1 ? v+1 : MAX_VERTEX_T)), v+1);
        }
    }
}

TEST(TopologyTreeSuite, path_query_test) {
    vertex_t n = 16;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    int id = 0;
    int d = 0;
    TopologyTree<int> tree(n, qt, f, id, d);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        // for (vertex_t v = 0; v <= i+1; v++) {
        //     ASSERT_EQ(tree.subtree_query(v, (v > 0 ? v-1 : MAX_VERTEX_T)), i+2-v);
        //     ASSERT_EQ(tree.subtree_query(v, (v <= i ? v+1 : MAX_VERTEX_T)), v+1);
        // }
    }
}
