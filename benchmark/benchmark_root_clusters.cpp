#include "benchmark.h"
#include "util.h"
#include "topology_tree.h"
#include "rc_tree.h"

void print_map(std::map<int,long long>  &m);
int total_root_clusters(std::map<int,long long>  &m);

int main(int argc, char** argv){
    std::vector<u_int32_t> n_list = {1000};
    std::tuple<std::string, std::function<std::vector<Update>(vertex_t, long)>, bool, int> test_cases[] = {
    {"Linked List", dynamic_tree_benchmark::linked_list_benchmark, false, 1},
    {"Binary Tree", dynamic_tree_benchmark::binary_tree_benchmark, false, 1},
    {"Random Degree 3", dynamic_tree_benchmark::random_degree3_benchmark, false, 1},
    };

  for (auto n : n_list) {
    // Update speed while supporting queries benchmark
    for (auto test_case : test_cases) {
      std::string test_case_name = std::get<0>(test_case);
      auto update_generator = std::get<1>(test_case);
      bool ternarize = std::get<2>(test_case);
      int num_trials = std::get<3>(test_case);
      std::vector<std::vector<Update>> update_sequences;
      for (int i = 0; i < num_trials; i++) {
        std::vector<Update> updates = update_generator(n, rand());
        // std::vector<Update> updates = dynamic_tree_benchmark::dense_tree_update_generator(update_generator(n, rand()));
        update_sequences.push_back(updates);
    }
    for(auto updates : update_sequences){
        TopologyTree<int, int> topT(n);
        RCTree<int> rcT(n);
        for (Update update : updates) {
            if (update.type == INSERT) {
                topT.link(update.edge.src, update.edge.dst);
                rcT.link(update.edge.src, update.edge.dst);
            } else if (update.type == DELETE) {
                topT.cut(update.edge.src, update.edge.dst);
                rcT.cut(update.edge.src, update.edge.dst);
            } else {
                std::cerr << "Invalid update type: " << update.type << std::endl;
                std::abort();
            }
        }
        std::cout << "Topology Tree Stats \n";
        std::cout << total_root_clusters(topology_root_clusters_histogram) << "\n";
        std::cout << "RC Tree Stats \n";
        std::cout << total_root_clusters(rc_root_clusters_histogram) << "\n";
    }
    }
  }
}

void print_map(std::map<int,long long>  &m) {
	for (auto it = m.cbegin(); it != m.cend(); it++) {	
		std::cout << (it->first) << ":"  << it->second  << "\n";
	}
}

int total_root_clusters(std::map<int,long long>  &m){
    int total = 0;
    for (auto it = m.cbegin(); it != m.cend(); it++) {	
        total += (it->first * it->second);
	}
    return total;
}