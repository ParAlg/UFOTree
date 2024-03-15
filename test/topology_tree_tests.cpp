#include <gtest/gtest.h>
#include "../include/topology_tree.h"


TEST(TopologyTreeSuite, example_test) {
    vertex_t num_vertices = 10;
    QueryType qt = PATH;
    auto f = [](int x, int y)->int{return x + y;};

    TopologyTree<int> tree(num_vertices, qt, f);
    tree.link(0,1);
    tree.cut(0,1);
}
