#include <gtest/gtest.h>
#include "../include/topology_tree.h"


TEST(TopologyTreeSuite, incremental_linkedlist_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v <= i+1; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
    }
}

TEST(TopologyTreeSuite, incremental_binarytree_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f);

    for (vertex_t i = 0; i < (n-1)/2; i++) {
        tree.link(i,2*i+1);
        tree.link(i,2*i+2);
        for (vertex_t u = 0; u < 2*i+2; u++) for (vertex_t v = u+1; v <= 2*i+2; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
    }
    if (n%2 == 0) tree.link((n-1)/2,n-1);
    for (vertex_t u = 0; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
}

TEST(TopologyTreeSuite, decremental_linkedlist_test) {
    vertex_t n = 256;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(n, qt, f);

    for (vertex_t i = 0; i < n-1; i++) {
        tree.link(i,i+1);
    }
    for (vertex_t i = 0; i < n-1; i++) {
        tree.cut(i,i+1);
        for (vertex_t u = 0; u < i+1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_FALSE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " connected.";
        for (vertex_t u = i+1; u < n-1; u++) for (vertex_t v = u+1; v < n; v++)
            ASSERT_TRUE(tree.connected(u,v)) << "Vertex " << u << " and " << v << " not connected.";
    }
}
