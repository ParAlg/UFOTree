#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <parlay>
#include <unordered_set>
#include<vector>

typedef uint32_t vertex_t;

#define GET_NEIGHBOR(source, cluster) cluster.boundary_nodes[0] != source ? cluster.boundary_nodes[0] : cluster.boundary_nodes[1]

template<typename aug_t>
struct RCCluster {
  aug_t aug_val;
  std::vector<int> boundary_nodes;
  RCCluster* parent;

  RCCluster::RCCluster(aug_t aug_val):aug_val(aug_val);

};

template<typename aug_t>
class RCTree {
  parlay::sequence<RCCluster<aug_t>> clusters;
  parlay::sequence<RCCluster<aug_t>> leaf_clusters;
  parlay::sequence<parlay::sequence<parlay::sequence<RCCluster>>> adj; // RCCluster pointer not int
  std::unordered_set<int> affected; 
  int degree_bound, n; //n is number of vertices in the tree.

public:
  RCTree::RCTree(int n, int degree_bound);
  /*RCTree::RCTree(vector<int[3]> tree, int n, int degree_bound); */

  RCTree::link(int u, int v, int weight);
  RCTree::cut(int u, int v);

};

// Constructor. Presently assume a degree bound of 3.
RCTree::RCTree(int n, int degree_bound) { 
  degree_bound = degree_bound;
  n = n;
  adj.resize(n);
}

/* ------- HELPER METHODS ------- */
int RCTree::get_degree(int v, int round) {

  int degree = 0;

  for(auto cluster : adj[round][v])
  if(cluster.boundaries.size() == 2) degree++;

  return degree;
}

bool RCTree::contracts(int vertex, int round) {
  // Take in vertex and round number to determine if the vertex contracts in that
  // round.
  if(round < adj.size() && adj[round + 1][v].size() == 0){
    return true;
  }
  return false;
}
/* ------------------------------- */
void RCTree::spread_affection(int round) {

  //Vertices that affection is spread to from orginal affected vertices
  std::unordered_set<int> new_affected; 
  auto induced_tree = adj[round];

  for(auto vertex : affected) {
    //Iterate through all adjacencies of affected vertex
    for(int i = 0; i < induced_tree[vertex].size(); i++) {

      cluster = induced_tree[vertex][i];

      if(get_degree(cluster, round) == 2) {
        int neighbor = GET_NEIGHBOR(vertex, cluster);
        if(!new_affected.contains(neighbor) && spread_by_dependence(neighbor, round, &new_affected)){
          new_affected.insert(neighbor);
        }
      } 
    }
  } 
} 
affected.merge(new_affected);
}

bool RCTree::spread_by_dependence(int vertex, int round, std::unordered_set<int> &current_affected){

  bool spreads = true;
  for(auto neighbor : induced_tree[vertex]){
    if(!contracts(neighbor) && (current_affected.contains(neighbor)) {
      spreads = false;
      break;
    }
  }
  return spreads;
}
std::unordered_set<int> RCTree::MIS(int round){

  //pick a random vertex
  //Include vertex in MIS
  //Delete Neighbors
  //Return
  std::unordered_set<int> maximal_set(affected);
  for(auto random_vertex : affected){

    // If the random vertex is not int the maximal_set then the vertex has already
    // been removed as a neighbor of some other vertex and we can't do any ops on it.
    if(!maximal_set.contains(random_vertex)) continue;
    // A few cases to deal with :
    // 1) If the vertex is next to a contracting unaffected vertex it is not to be considered.
    // 2) If the vertex cannot contract, it must remain in the set of affected vertices.
    
    // Case if next to unaffected contracting vertex.
    boolean adjacent_unaffected_contracting = false;
    for(auto cluster : adj[round][random_vertex]){

      if(get_degree(cluster) == 2){
        int neighbor = GET_NEIGHBOR(random_vertex, cluster);
        if(!maximal_set.contains(neighbor) && contracts(neighbor)){
          maximal_set.remove(random_vertex); // Remove from MIS, if so.
          adjacent_unaffected_contracting = true;
        }
      }
    }

    // If not case 1 then this vertex is part of indep. set. Remove all neighbors.
    if(!adjacent_unaffected_contracting){
      for(auto cluster: adj[round][random_vertex]){
        int neighbor = GET_NEIGHBOR(random_vertex, cluster);
        if(maximal_set.contains(neighbor)) maximal_set.remove(neighbor);
      }
    }
  }
  
  return maximal_set;
}

RCTree::rake(int vertex, int round){

}

RCTree::compress(int vertex, int round){

}

RCTree::link(int u, int v, int weight) {
  // Add an edge between 2 Trees of the RC forest.
  // Assume that the vertices u and v always exist.
  // and the edge being added does not yet violate
  // any tree properties.

  parlay::sequence<parlay::sequence<RCCluster>> base_tree = adj[0];

  RCCluster new_edge(weight);
  new_edge.boundary_vertexes[0] = u;
  new_edge.boundary_vertexes[1] = v;

  base_tree[u].insert(new_edge); //v can be > t. Make it ternary.
  base_tree[v].insert(new_edge);

  affected.insert(u);          // Insert initial affected vertices.
  affected.insert(v);          

  update(); 
}

RCTree::cut(int u, int v) {
  //Delete an edge between 2 vertices of an RC Tree.
  //If no edge (u,v) exists, throw an exception.
}

RCTree::update() {
  // Based on set of affected vertices, find an MIS of the affected vertices, which induce
  // a subtree on the original tree and recontract them accordingly. Then determine new
  // set of affected vertices and recurse.
  
  // Spread affection and redo contractions till there are affected
  // vertices.
  int round = 0
  while(!affected.empty()){
    spread_affection(round);
    std::unordered_set<vertex_t> maximal_set = MIS(round);
    for(auto vertex : maximal_set){
      //Redo contractions for all vertices in maximal set.
      //Remove vertices from set of affected vertices.
    }
  }
}
