#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/top_tree.h"

TEST(Top_tree_suite, constructor_test){
    TopTree<int> t(3, PATH, [] (int a, int b){return std::min(a,b);}, std::numeric_limits<int>::max(), std::numeric_limits<int>::max());
    ASSERT_EQ(t.t.num_vertices, 3);
}