#include <gtest/gtest.h>
#include "graph_utils.h"


int main(int argc, char** argv) {
  std::string file_names[] = {
      "/ssd1/zhongqi/graphdata/sym/com-orkut_sym.bin"
  };

  for (auto file_name : file_names) {
      auto G = graph_utils::break_sym_graph_from_bin(file_name);
      graph_utils::print_graph_stats(G);
      std::cout << "diameter = " << graph_utils::get_graph_diameter(G);
  }

  testing::InitGoogleTest(&argc, argv);
  int ret = RUN_ALL_TESTS();
  return ret;
}
