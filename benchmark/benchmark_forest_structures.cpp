#include "forest_strucure.h"
#include <parlay/sequence.h>


using namespace parlay;

int main(int argc, char** argv) {
    std::cout << "Benchmarking forest structures." << std::endl;
    SequentialForest F;
    sequence<vertex_t> vertices;
    F.insert_vertices(vertices);
    F.delete_vertices(vertices);
    sequence<Edge> edges;
    F.insert_edges(edges);
    F.delete_edges(edges);
    F.check_edges(edges);
}
