#include "forest_strucure.h"
#include "sequential_forest.h"
#include <parlay/sequence.h>


using namespace parlay;

void benchmark_forest_structure(ForestStructure* F);

int main(int argc, char** argv) {
    std::cout << "Benchmarking forest structures." << std::endl;
    ForestStructure* F;

    SequentialForest SF;
    F = &SF;
    benchmark_forest_structure(F);
}

void benchmark_forest_structure(ForestStructure* F) {
    sequence<vertex_t> vertices;
    F->insert_vertices(vertices);
    F->delete_vertices(vertices);

    sequence<Edge> edges;
    F->insert_edges(edges);
    F->delete_edges(edges);
    F->check_edges(edges);

    F->get_degree(0);

    sequence<vertex_t> parents = F->get_parents(vertices);
    sequence<pair<vertex_t,vertex_t>> parent_counts = F->count_parents(vertices);
    F->set_parents(vertices, vertices);
    F->unset_parents(vertices);

    sequence<vertex_t> child_counts;
    F->add_children(child_counts);
    F->get_child_count(0);
}
