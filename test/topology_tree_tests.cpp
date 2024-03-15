#include <gtest/gtest.h>
#include "../include/topology_tree.h"


TEST(TopologyTreeSuite, example_test) {
    auto f = [](int x, int y)->int{return x + y;};
    TopologyTree<int> tree(10, PATH, f);
}
