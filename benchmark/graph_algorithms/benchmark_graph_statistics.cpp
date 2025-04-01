#include "graph_benchmark.h"


int main(int argc, char** argv) {
    std::string file_names[] = {
        "/ssd1/zhongqi/graphdata/sym/com-youtube_sym.bin",
        "/ssd1/zhongqi/graphdata/sym/com-orkut_sym.bin",
        "/ssd1/zhongqi/graphdata/sym/soc-LiveJournal1_sym.bin",
        "/ssd1/zhongqi/graphdata/sym/twitter_sym.bin",
        "/ssd1/zhongqi/graphdata/sym/friendster_sym.bin",
    };

    for (auto file_name : file_names) {
        auto G = graph_utils::break_sym_graph_from_bin(file_name);
        auto BFS_F = graph_utils::BFS_forest(G);
        auto RIS_F = graph_utils::RIS_forest(G);

        std::cout << "[ GRAPH: " << file_name.substr(file_name.find_last_of('/')+1) << " ]" << std::endl;
        graph_utils::print_graph_stats(G);
        std::cout << "Graph CC Count = " << graph_utils::get_component_count(G) << std::endl;
        std::cout << "Graph Diameter = " << graph_utils::get_graph_diameter(G) << std::endl;
        std::cout << "BFS Forest Diameter = " << graph_utils::get_forest_diameter(BFS_F, G.size()) << std::endl;
        std::cout << "RIS Forest Diameter = " << graph_utils::get_forest_diameter(RIS_F, G.size()) << std::endl;
    }
}