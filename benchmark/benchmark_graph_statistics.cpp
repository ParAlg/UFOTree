#include "graph_benchmark.h"


int main(int argc, char** argv) {
    std::string file_names[] = {
        "/ssd1/zhongqi/graphdata/sym/com-orkut_sym.bin"
    };

    for (auto file_name : file_names) {
        auto G = graph_utils::break_sym_graph_from_bin(file_name);
        auto F = graph_utils::BFS_forest(G);

        graph_utils::print_graph_stats(G);
        std::cout << "BFS Forest Diameter = " << graph_utils::get_forest_diameter(F, G.size()) << std::endl;
        std::cout << "Graph Diameter      = " << graph_utils::get_graph_diameter(G) << std::endl;
    }
}