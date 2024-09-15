#include <gtest/gtest.h>


int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
  testing::GTEST_FLAG(filter) = "ParallelUFO*";
  int ret = RUN_ALL_TESTS(); 
  return ret; 
}
