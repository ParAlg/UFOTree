#include <cstdio>
#include <cstdlib>
#include <iterator>
#include <parlay/sequence.h>
#include <unordered_set>
#include<vector>
#include<parlay/monoid.h>

typedef uint32_t vertex_t;

#define GET_NEIGHBOR(source, cluster) cluster.boundary_nodes[0] != source ? cluster.boundary_nodes[0] : cluster.boundary_nodes[1]

template<typename aug_t>
struct RCCluster {
public:
  aug_t aug_val;
  vertex_t representative_vertex = -1; //Set to -1 to recognize uncontracted binary clusters.
  std::vector<int> boundary_vertexes;
  RCCluster* parent;
    RCCluster(aug_t _aug_val){ aug_val = _aug_val;}
    
    RCCluster(aug_t _aug_val, vertex_t _representative_vertex){
      aug_val = _aug_val;
      representative_vertex = _representative_vertex;
    }
};

template<typename aug_t>
class RCTree {
public: 
  int degree_bound, n; //n is number of vertices in the tree.
  parlay::sequence<RCCluster<aug_t>> clusters;
  parlay::sequence<RCCluster<aug_t>> leaf_clusters;
  std::unordered_set<int> affected; 
  parlay::sequence<RCCluster<aug_t>***> adj; // RCCluster pointer not int - Made public for testing, will actually be private.
  void add_neighbor(RCCluster<aug_t>** adjacencies, RCCluster<aug_t>* cluster);
  int get_degree(int v, int round);
  bool contracts(vertex_t v, int round);
  void spread_affection(int round);
  bool spread_by_dependence(int vertex, int round, std::unordered_set<int> &current_affected);
  std::unordered_set<int> MIS(int round);
  bool contains_rep(vertex_t neighbor, vertex_t v, int round);
  void rake(vertex_t vertex, int round);
  void compress(vertex_t vertex, int round);
  void update();
  int neighbor_count(RCCluster<aug_t>** neighbors);
  RCTree(int n, int degree_bound);
  /*RCTree::RCTree(vector<int[3]> tree, int n, int degree_bound); */
  void link(int u, int v, int weight);
  void cut(int u, int v);
};

template<typename aug_t>
// Constructor. Presently assume a degree bound of 3.
RCTree<aug_t>::RCTree(int _n, int _degree_bound) { 
  degree_bound = _degree_bound;
  n = _n;
  adj.push_back(new RCCluster<aug_t>** [_n]);

  for(int i = 0; i < n; i++){ adj[0][i] = (RCCluster<int>**) calloc(degree_bound, sizeof(RCCluster<int>*));}
}

/* ------- HELPER METHODS ------- */
template<typename aug_t>
int RCTree<aug_t>::get_degree(int v, int round) {

  int degree = 0;

  for(int i = 0; i < degree_bound; i++){
    auto cluster = adj[round][v][i];
    if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) degree++;
  }

  return degree;
}

template<typename aug_t>
int RCTree<aug_t>::neighbor_count(RCCluster<aug_t>** neighbors){
  int count = 0;
  for(int i = 0; i < degree_bound; i++){
    auto neighbor_ptr = neighbors[i];
    if(neighbor_ptr != nullptr){
      count += 1;
    }
  }
  return count;
}

template<typename aug_t>
bool RCTree<aug_t>::contracts(vertex_t v, int round) {
  // Take in vertex and round number to determine if the vertex contracts in that
  // round.
  if(round < adj.size() - 1 && neighbor_count(adj[round + 1][v]) == 0){
    return true;
  }
  return false;
}

template<typename aug_t>
void RCTree<aug_t>::add_neighbor(RCCluster<aug_t>** adjacencies, RCCluster<aug_t> *cluster){
  for(int i = 0; i < degree_bound; i++){
    if(adjacencies[i] == nullptr){
      adjacencies[i] = cluster;
      break;
    }
  }
}
/* ------------------------------- */
template<typename aug_t>
void RCTree<aug_t>::spread_affection(int round) {
  /* Determine the vertices we spread affection to via dependence as a result of changes
   * to the input and add them to the affected set, prior to building an MIS of the
   * affected vertices.
   */

  //Vertices that affection is spread to from orginal affected vertices
  std::unordered_set<int> new_affected(affected); 
  auto induced_tree = adj[round];

  for(auto vertex : affected) {
    //Iterate through all adjacencies of affected vertex
    for(int i = 0; i < degree_bound; i++) {

      RCCluster<aug_t> cluster = *induced_tree[vertex][i];

      if(get_degree(cluster, round) == 2) {
        int neighbor = GET_NEIGHBOR(vertex, cluster); // Neighbor of the affected vertex.
        if(!new_affected.count(neighbor) && spread_by_dependence(neighbor, round, &new_affected)){
          new_affected.insert(neighbor);
        }
      } 
    }
  } 
}

template<typename aug_t>
bool RCTree<aug_t>::spread_by_dependence(int vertex, int round, std::unordered_set<int> &current_affected){
  /* Determine if a vertex becomes affected as result of contracting neighbors added to the affected set.
   * Return true if the vertex is now affected, false otherwise.*/ 
  auto induced_tree = adj[round];
  bool spreads = true;
  for(auto neighbor : induced_tree[vertex]){ //Neighbors of the neighbor of the affected vertex.
    if(!contracts(*neighbor) && (current_affected.count(*neighbor))) {
      spreads = false;
      break;
    }
  }
  return spreads;
}

template<typename aug_t>
std::unordered_set<int> RCTree<aug_t>::MIS(int round){

  //pick a random vertex
  //Include vertex in MIS
  //Delete Neighbors
  //Return
  std::unordered_set<int> maximal_set(affected);
  for(auto random_vertex : affected){

    // If the random vertex is not int the maximal_set then the vertex has already
    // been removed as a neighbor of some other vertex and we can't do any ops on it.
    if(!maximal_set.count(random_vertex)) continue;
    // A few cases to deal with :
    // 1) If the vertex is next to a contracting unaffected vertex it is not to be considered.
    // 2) If the vertex cannot contract, it must remain in the set of affected vertices.
    
    // Case if next to unaffected contracting vertex.
    bool adjacent_unaffected_contracting = false;
    for(auto cluster : adj[round][random_vertex]){

      if(get_degree(cluster) == 2){
        int neighbor = GET_NEIGHBOR(random_vertex, cluster);
        if(!maximal_set.count(neighbor) && contracts(neighbor)){
          maximal_set.erase(random_vertex); // Remove from MIS, if so.
          adjacent_unaffected_contracting = true;
          break;
        }
      }
    }

    // If not case 1 and vertex degree < 3, then this vertex is part of indep. set. 
    // Remove all neighbors.
    if(!adjacent_unaffected_contracting && get_degree(random_vertex) != 3){
      for(auto cluster: adj[round][random_vertex]){
        int neighbor = GET_NEIGHBOR(random_vertex, cluster);
        if(maximal_set.count(neighbor)) maximal_set.erase(neighbor);
      }
    }
  }
  
  return maximal_set;
}
template<typename aug_t>
void RCTree<aug_t>::link(int u, int v, int weight) {
  // Add an edge between 2 Trees of the RC forest.
  // Assume that the vertices u and v always exist.
  // and the edge being added does not yet violate
  // any tree properties.

  auto base_tree = adj[0];

  RCCluster new_edge(weight);
  new_edge.boundary_vertexes[0] = u;
  new_edge.boundary_vertexes[1] = v;

  add_neighbor(base_tree[u], &new_edge); //v can be > t. Make it ternary.
  add_neighbor(base_tree[v], &new_edge);

  affected.insert(u);          // Insert initial affected vertices.
  affected.insert(v);          

  update(); 
}

template<typename aug_t>
void RCTree<aug_t>::cut(int u, int v) {
  //Delete an edge between 2 vertices of an RC Tree.
  //If no edge (u,v) exists, throw an exception.
}

template<typename aug_t>
bool RCTree<aug_t>::contains_rep(vertex_t neighbor, vertex_t v, int round){
  auto adjacencies = adj[round][neighbor];
  RCCluster<aug_t> to_return;
  for(auto neighboring_cluster : adjacencies){
    if(neighboring_cluster.representative_vertex == -1 && 
        neighboring_cluster.boundary_nodes.count(v)){
      to_return = neighboring_cluster;
      break;
    }

    if(neighboring_cluster.representative_vertex == v){
      to_return = neighboring_cluster;
      break;
    }
  }
  return to_return;
}

template<typename aug_t>
void RCTree<aug_t>::rake(vertex_t vertex, int round){
  //Get boundary vertex
  auto neighbors = adj[round][vertex];
  vertex_t neighbor; //To be initialized later, helpful for doing contractions.
  RCCluster<aug_t> new_cluster(0, vertex);  //New augmented cluster to be used for replacement. - Monoid struct in parlay include.

  // Go through all neighboring clusters and aggregate values onto self.
  for(auto neighboring_cluster : neighbors) {
    new_cluster.aug_val += neighboring_cluster.aug_val; // Needs to be a associate binary op passed in by the user.
    if(neighboring_cluster.boundary_vertexes.size() == 2){
      neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
    }
  }

  new_cluster.boundary_vertexes.insert(neighbor);
  //Check neighbor's clusters to see if we exist as a rep vertex and if so
  //compare current cluster to it.
  RCCluster<aug_t> representative_cluster = contains_rep(neighbor, vertex, round + 1);

  // Add other vertex I am raking onto into affected set if 
  // the cluster of the vertex being raked in the next round
  // is not a unary cluster. 
  // QUESTION : DO I NEED TO ADD IT IF VALUES NOT EQUIVALENT AS WELL - YES, to recontract with new values.
  affected.insert(neighbor);

  adj[round + 1][neighbor].remove(representative_cluster);
  adj[round + 1][neighbor].insert(new_cluster);

  //Now that the vertex has contracted in this round, we must clear everything
  //else that happened when this vertex was uncontracted in subsequent rounds.
  //QUESTION : Can this cause issues with other contractions that were adjacent
  //to this one in subsequent rounds trying to contract themselves onto this.
}

template<typename aug_t>
void RCTree<aug_t>::compress(vertex_t vertex, int round){
  
}

template<typename aug_t>
void RCTree<aug_t>::update() {
  // Based on set of affected vertices,.count an MIS of the affected vertices, which induce
  // a subtree on the original tree and recontract them accordingly. Then determine new
  // set of affected vertices and recurse.
  
  // Spread affection and redo contractions till there are affected
  // vertices.
  int round = 0;
  while(!affected.empty()){
    spread_affection(round);
    std::unordered_set<vertex_t> maximal_set = MIS(round);
    for(auto vertex : maximal_set){
      //Redo contractions for all vertices in maximal set.
      //Remove vertices from set of affected vertices.
      if(get_degree(vertex, round) == 1) {
        rake(vertex, round);
      } else {
        compress(vertex, round);
      }
    }
  }
}
