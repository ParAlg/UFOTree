#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <stdexcept>
#include <unordered_set>
#include "../include/spaa_rc_tree.h"

TEST(ParallelRCTreeSuite, test_constructor){
  parRCTree<int> t(3);
  /*for(int i = 0; i < t.clusters.size(); i++){
    t.clusters[i].print();
  }*/
  ASSERT_EQ(t.clusters.size(), 3);
}

TEST(ParallelRCTreeSuite, test_batch_insert){

}
