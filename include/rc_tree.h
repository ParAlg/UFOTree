#pragma once
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <parlay/sequence.h>
#include <stdexcept>
#include <sys/param.h>
#include <unordered_set>
#include<vector>
#include "../util/types.h"
#include "../util/util.h"
#include "ternarizable_interface.h"

//#define COLLECT_ROOT_CLUSTER_STATS

#define GET_NEIGHBOR(source, cluster) cluster.boundary_vertexes[0] != source ? cluster.boundary_vertexes[0] : cluster.boundary_vertexes[1]

#define DEGREE_BOUND 3 // CHANGE THIS IN CASE YOU WANT A DIFFERENT DEGREE BOUND FOR YOUR PURPOSES.

#ifdef COLLECT_ROOT_CLUSTER_STATS
    static std::map<int, int> rc_root_clusters_histogram;
#endif


namespace ufo {

template<typename aug_t>
struct RCCluster {
public:
  aug_t aug_val;
  vertex_t boundary_vertexes[2];

  RCCluster<aug_t>* parent;
  vertex_t rep_vertex; // For testing purposes mainly. Maybe necessary for subtree queries?

  RCCluster(aug_t _aug_val){ aug_val = _aug_val; parent = nullptr; rep_vertex = MAX_VERTEX_T;boundary_vertexes[0] = MAX_VERTEX_T; boundary_vertexes[1] = MAX_VERTEX_T;}
  RCCluster(aug_t _aug_val, RCCluster<aug_t>* _parent){
    aug_val = _aug_val;
    parent = _parent;
    rep_vertex = MAX_VERTEX_T;
    boundary_vertexes[0] = MAX_VERTEX_T; boundary_vertexes[1] = MAX_VERTEX_T;
  }
  int bv_size();
  void add_boundary(vertex_t v);
  
  bool operator==(const RCCluster<aug_t> &other){
    return (aug_val == other.aug_val && 
            boundaries_equal(this, &other) &&
            parent == other.parent && rep_vertex == other.rep_vertex);
  }
};

template<typename aug_t>
class RCTree : ITernarizable<aug_t> {
public:
  // Data structure fields.
  vertex_t n; //n is number of vertices in the tree.
  aug_t id, d_val; //identity element, default value
  std::function<aug_t(aug_t, aug_t)> f; // Associative function
  std::vector<vertex_t> roots; //for debugging purposes only.
  std::vector<vertex_t> affected;
  std::vector<std::vector<RCCluster<aug_t>**>> contraction_tree; 
  std::vector<RCCluster<aug_t>*> representative_clusters;   
  std::vector<RCCluster<aug_t>*> to_delete;
  std::vector<vertex_t> maximal_set;
  std::vector<vertex_t> new_affected;
  std::vector<vertex_t> uncontracting;
  std::vector<int> round_contracted; 

  // Functions
  void is_valid_MIS(int round);
  void add_neighbor(int round, RCCluster<aug_t>* cluster, vertex_t v);
  int get_degree(vertex_t v, int round);
  bool contracts(vertex_t v, int round, bool version);
  void spread_affection(int round);
  bool affected_by_dependence(vertex_t vertex, int round, std::vector<vertex_t> *current_affected);
  void MIS(int round);
  void rake(vertex_t vertex, int round);
  void compress(vertex_t vertex, int round);
  void update();
  int neighbor_count(RCCluster<aug_t>** neighbors);
  void finalize(vertex_t vertex, int round);
  RCCluster<aug_t>* get_root(RCCluster<aug_t>* cluster);
  void clear_deleted_clusters(vertex_t vertex, int round);
  void is_valid_induced_tree(int round);
  bool connected(vertex_t u, vertex_t v);
  bool is_edge(RCCluster<aug_t> *cluster);
  bool edge_exists(vertex_t u, vertex_t v);
  RCTree<aug_t>** get_neighbors(vertex_t v);
  size_t space(); // For space benchmarks.

  // Helper functions for "vector-set"
  bool vector_contains(std::vector<vertex_t>* vec, vertex_t v);
  void vec_erase(std::vector<vertex_t>* vec, vertex_t v);
  void vec_insert(std::vector<vertex_t>* vec, vertex_t v);
  
  bool vector_contains(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v);
  void vec_erase(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v);
  void vec_insert(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v);
  // Queries
  aug_t path_query(vertex_t u, vertex_t v);

  // Overrides for ternarizable_interface methods
  short get_degree(vertex_t v) override {return get_degree(v, 0);}
  std::pair<vertex_t, aug_t> retrieve_v_to_del(vertex_t v) override {
    return std::pair(GET_NEIGHBOR(v, (*contraction_tree[v][0][0])), contraction_tree[v][0][0]->aug_val);
  }

  /*RCTree::RCTree(vector<int[3]> tree, int n, int DEGREE_BOUND); */
  // These are only ones that should really be public.
  RCTree(int _n, QueryType q = PATH, 
         std::function<aug_t(aug_t, aug_t)> f = [] (aug_t x, aug_t y) {return x + y;}, 
         aug_t id = 0, aug_t d_val = 0);

  RCTree(int _n, QueryType _q, std::function<aug_t(aug_t, aug_t)> f_v, std::function<aug_t(aug_t, aug_t)> f_e, aug_t id_v, aug_t id_e, aug_t d_val_v, aug_t d_val_e);
  ~RCTree();
  void link(vertex_t u, vertex_t v, aug_t weight = 1);
  void cut(vertex_t u, vertex_t v);
};

template<typename aug_t>
// Constructor. 
RCTree<aug_t>::RCTree(int _n, QueryType _q, 
                      std::function<aug_t(aug_t, aug_t)> _f,
                      aug_t _id, aug_t _d_val) { 

  n = _n; f = _f; id = _id; d_val = _d_val; 
  // Initialize round 0 adjacency list for each vertex in the tree.
  for(int i = 0; i < _n; i++){
    contraction_tree.push_back(std::vector<RCCluster<aug_t>**>());
  }
  representative_clusters.resize(n); 
  round_contracted.resize(n); 
  for(int i = 0; i < n; i++){ 
    contraction_tree[i].push_back(new RCCluster<aug_t>*[DEGREE_BOUND]{nullptr});
    representative_clusters[i] = (RCCluster<aug_t>*) new RCCluster<aug_t>(id);
    representative_clusters[i]->rep_vertex = i;
    round_contracted[i] = 0;
  }
}

// FOR USE IN Ternarized_Tree CLASS ONLY
template<typename aug_t>
RCTree<aug_t>::RCTree(int _n, QueryType _q, std::function<aug_t(aug_t, aug_t)> f_v, std::function<aug_t(aug_t, aug_t)> f_e, aug_t id_v, aug_t id_e, aug_t d_val_v, aug_t d_val_e)
  : RCTree<aug_t>(_n, _q, f_e, id_e, d_val_e) {}
// Destructor
template<typename aug_t>
RCTree<aug_t>::~RCTree(){
  for(int i = 0; i < representative_clusters.size(); i++){
    delete representative_clusters[i];
  }
  for(int i = 0; i < contraction_tree.size(); i++){
    for(int j = 0; j < contraction_tree[i].size(); j++){
      delete[] contraction_tree[i][j];
    }
  }
}
/* ------- HELPER METHODS ------- */
template<typename aug_t>
size_t RCTree<aug_t>::space(){
  size_t memory = sizeof(RCTree<aug_t>);
  for(auto ptr : representative_clusters){
    memory += (sizeof(*ptr) + sizeof(ptr));
  }
  for(int i = 0; i < contraction_tree.size(); i++){
    for(int j = 0; j < contraction_tree[i].size(); j++){
      for(int k = 0; k < DEGREE_BOUND; k++){
        memory += sizeof(contraction_tree[i][j][k]);
      }
    }
    memory += (contraction_tree[i].size() * sizeof(RCCluster<aug_t>*));
  }
  memory += (contraction_tree.size() * sizeof(std::vector<RCCluster<aug_t>**>));
  memory += round_contracted.size() * sizeof(int);
  return memory;
}

template<typename aug_t>
int RCTree<aug_t>::get_degree(vertex_t v, int round) {
  // Get the degree of a vertex at a particular round i.e. no. of binary clusters in the
  // vertex's adjacency list at a round.
  int degree = 0;
  for(int i = 0; i < DEGREE_BOUND; i++){
    auto cluster = contraction_tree[v][round][i];
    if(is_edge(cluster)) degree++;
  }
  return degree;
}

template<typename aug_t>
bool RCTree<aug_t>::contracts(vertex_t v, int round, bool version) {
  // Parameters are the vertex, round number and which variant
  // of the method to run. Version 0 (represented by bool version variables) will
  // return true if vertex contracts exactly in that round. Version 1 will return
  // true if vertex contracts in that round or any previous round.
  if(version == 0){
    return round_contracted[v] == round; 
  }
  return round_contracted[v] <= round;
}

template<typename aug_t>
void RCTree<aug_t>::add_neighbor(int round, RCCluster<aug_t>* cluster, vertex_t contracted_v){
  /* Add a cluster into the adjacency list of its boundary vertex(es) at a particular
   * round*/
  for(int i = 0; i < cluster->bv_size();i++){
    //Get adjacency list of a boundary vertex of the cluster at a round.
    auto adjacencies = contraction_tree[cluster->boundary_vertexes[i]][round]; 
    auto vertex_added = false;

    for(int i = 0; i < DEGREE_BOUND; i++){
      if(adjacencies[i] == nullptr){
        adjacencies[i] = cluster;
        vertex_added = true;
        break;
      }
    }
    /* This method is only called when it is possible to add a vertex into the adjacency list of its
     * boundary vertices, if this does not happen then there is an error.*/
    if(!vertex_added){
      throw std::invalid_argument("Vertex could not be added to the tree");
    }
  }
}


// Method to check if boundary vertices of 2 clusters are equal.
template<typename aug_t>
bool boundaries_equal(RCCluster<aug_t>* neighbor, RCCluster<aug_t>* cluster){
  return (neighbor->boundary_vertexes[0] == cluster->boundary_vertexes[0] &&
          neighbor->boundary_vertexes[1] == cluster->boundary_vertexes[1]) ||
         (neighbor->boundary_vertexes[0] == cluster->boundary_vertexes[1] &&
          neighbor->boundary_vertexes[1] == cluster->boundary_vertexes[0]);
}

template<typename aug_t>
RCCluster<aug_t>* RCTree<aug_t>::get_root(RCCluster<aug_t>* cluster){
  // Returns the root of a cluster in the RC Tree if it exists.
  auto curr = cluster;
  while(curr->parent != nullptr){
    curr = curr->parent; 
  }
  return curr;
}

template<typename aug_t>
bool RCTree<aug_t>::connected(vertex_t u, vertex_t v){
  // Checks to see if 2 vertices are connected i.e. have a path between
  // them by verifying if they are in the same component in the RC Tree.
  return get_root(representative_clusters[u]) == get_root(representative_clusters[v]);
}

template<typename aug_t>
bool RCTree<aug_t>::is_edge(RCCluster<aug_t> *cluster){
  // Check to see if an RCCluster is a binary cluster.
  return cluster != nullptr && cluster->bv_size() == 2;
}

template<typename aug_t>
bool RCTree<aug_t>::edge_exists(vertex_t u, vertex_t v){
  auto u_adjacencies = contraction_tree[u][0];
  for(int i = 0; i < DEGREE_BOUND; i++){
    if(is_edge(u_adjacencies[i])){
      if((GET_NEIGHBOR(u, (*u_adjacencies[i]))) == v){
        bool contains_cluster = false;
        for(int j = 0; j < DEGREE_BOUND; j++){
          if(contraction_tree[v][0][j] == u_adjacencies[i]) contains_cluster = true;
        }
        if(!contains_cluster){
          throw std::invalid_argument("Unidirectional edge found, edge_exists method.");
        }
        
        return true;
      }
    } 
  }
  return false;
}

template <typename aug_t>
void RCTree<aug_t>::clear_deleted_clusters(vertex_t vertex, int round){
  // This method is used to clear out clusters from the adjacency
  // lists of unaffected vertices that neighbor affected vertices
  // Clusters that represent an edge to an affected vertex or are
  // formed from the contraction of a vertex that is now affected
  // are replaced with nullptrs as they will be recomputed in the
  // present round.
  for(int i = 0; i < DEGREE_BOUND; i++){
    if(contraction_tree[vertex][round][i] != nullptr && 
        vector_contains(&to_delete, contraction_tree[vertex][round][i])){
      contraction_tree[vertex][round][i] = nullptr;
    }

    if(is_edge(contraction_tree[vertex][round][i])){
      auto neighbor = GET_NEIGHBOR(vertex, (*contraction_tree[vertex][round][i]));
      if(vector_contains(&affected, neighbor)){
        contraction_tree[vertex][round][i] = nullptr;
      }
    }
  }
}
template<typename aug_t>
RCTree<aug_t>** RCTree<aug_t>::get_neighbors(vertex_t v){
  return contraction_tree[v][0];
}

template<typename aug_t>
bool RCTree<aug_t>::vector_contains(std::vector<vertex_t>* vec, vertex_t v){
  return std::find(vec->begin(), vec->end(), v) != vec->end();
}

template<typename aug_t>
void RCTree<aug_t>::vec_erase(std::vector<vertex_t>* vec, vertex_t v){
  assert(vector_contains(vec, v));
  vec->erase(std::find(vec->begin(), vec->end(), v));
}

template<typename aug_t>
void RCTree<aug_t>::vec_insert(std::vector<vertex_t>* vec, vertex_t v){
  if(!vector_contains(vec, v)) vec->push_back(v);
}

template<typename aug_t>
bool RCTree<aug_t>::vector_contains(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v){
  return std::find(vec->begin(), vec->end(), v) != vec->end();
}

template<typename aug_t>
void RCTree<aug_t>::vec_erase(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v){
  assert(vector_contains(vec, v));
  vec->erase(std::find(vec->begin(), vec->end(), v));
}

template<typename aug_t>
void RCTree<aug_t>::vec_insert(std::vector<RCCluster<aug_t>*>* vec, RCCluster<aug_t>* v){
  if(!vector_contains(vec, v)) vec->push_back(v);
}
/* ------------------------------- */

template<typename aug_t>
void RCTree<aug_t>::spread_affection(int round) {
  /* Determine the vertices we spread affection to via dependence as a result of changes
   * to the input and add them to the affected set, prior to building an MIS of the
   * affected vertices.
   */

  //Vertices that affection is spread to from orginal affected vertices
  for(vertex_t v : affected) new_affected.push_back(v); 

  for(auto vertex : affected) {
    //Iterate through all adjacencies of affected vertex
    for(int i = 0; i < DEGREE_BOUND; i++) {

      RCCluster<aug_t>* neighboring_cluster = contraction_tree[vertex][round][i];

      if(is_edge(neighboring_cluster)) {
        vertex_t neighbor = GET_NEIGHBOR(vertex, (*neighboring_cluster)); // Neighbor of the affected vertex.
        //Check to see if neighbor is already affected and if all contracting neighbors are affected.
        if(!vector_contains(&new_affected, neighbor) && affected_by_dependence(neighbor, round, &new_affected)){ 
          vec_insert(&new_affected, neighbor); 
        }
      } 
    }
  }
  for(vertex_t v : new_affected) vec_insert(&affected, v); 
  new_affected.clear();
}

template<typename aug_t>
bool RCTree<aug_t>::affected_by_dependence(vertex_t vertex, int round, std::vector<vertex_t> *current_affected){
  /* Determine if a vertex becomes affected as result of contracting neighbors added to the affected set.
   * Return true if the vertex is now affected, false otherwise.*/ 
  if(contracts(vertex, round, 0)) return false;

  bool spreads = true;
  for(int i = 0; i < DEGREE_BOUND; i++){ //Neighbors of the neighbor of the affected vertex.
    auto neighbor_cluster = contraction_tree[vertex][round][i]; 
    if(is_edge(neighbor_cluster)){
      int neighbor = GET_NEIGHBOR(vertex, (*neighbor_cluster)); // Neighbor of the affected vertex.
      // If my neighbor does not contract
      if(contracts(neighbor, round, 0) && !vector_contains(current_affected, neighbor)) {
        spreads = false;
        break;
      }
    }
  }
  return spreads;
}

template<typename aug_t>
void RCTree<aug_t>::MIS(int round){
  /* Construct an MIS on the set of affected vertices at a particular round as follows: */

  // Initially MIS contains all vertices in affected set.
  // Pick a random vertex
  // Include vertex in MIS
  // Delete Neighbors
  // Recurse

  for(auto v : affected) vec_insert(&maximal_set, v); 

  for(auto random_vertex : affected){
    // If the random vertex is not in the maximal_set then the vertex has already
    // been removed as a neighbor of some other vertex and we can't do any ops on it.
    if(!vector_contains(&maximal_set, random_vertex)) continue;
    // A few cases to deal with :
    // 1) If the vertex is next to a contracting unaffected vertex it is not to be considered.
    // 2) If the vertex cannot contract, it must remain in the set of affected vertices.

    // Case if next to unaffected contracting vertex.
    bool adjacent_unaffected_contracting = false;
    for(int i = 0; i < DEGREE_BOUND; i++){
      auto cluster = contraction_tree[random_vertex][round][i];
      if(is_edge(cluster)) {
        int neighbor = GET_NEIGHBOR(random_vertex, (*cluster)); 

        //Neighbor not in affected set and neighbor contracts in this round.
        if(contracts(neighbor, round, 0) && !vector_contains(&affected, neighbor)){
          vec_erase(&maximal_set, random_vertex); // Remove from MIS, if so.
          adjacent_unaffected_contracting = true;
          break;
        }
      }
    }

    // If not case 1 and vertex degree < 3, then this vertex is part of indep. set. 
    // Remove all neighbors.
    if(!adjacent_unaffected_contracting && get_degree(random_vertex, round) != 3){
      for(int i = 0; i < DEGREE_BOUND; i++){
        auto cluster = contraction_tree[random_vertex][round][i];
        if(is_edge(cluster)){
          int neighbor = GET_NEIGHBOR(random_vertex, (*cluster));
          if(vector_contains(&maximal_set, neighbor)) vec_erase(&maximal_set, neighbor); 
        }
      }
    } else {
      if(vector_contains(&maximal_set, random_vertex)) vec_erase(&maximal_set, random_vertex);
    }
  }
}

template<typename aug_t>
void RCTree<aug_t>::link(vertex_t u, vertex_t v, aug_t weight) {
  // Add an edge between 2 Trees of the RC forest.
  assert(u >= 0 && v < n && !connected(u,v));

  // Create RC cluster to represent edge at round 0 between u and v.
  RCCluster<aug_t>* new_edge = (RCCluster<aug_t>*) new RCCluster<aug_t>(weight);
  new_edge->add_boundary(u);
  new_edge->add_boundary(v);

  add_neighbor(0, new_edge, MAX_VERTEX_T);
  vec_insert(&affected, u); 
  vec_insert(&affected, v); 
  roots.clear();
  update();
}

template<typename aug_t>
void RCTree<aug_t>::cut(vertex_t u, vertex_t v) {
  //Delete an edge between 2 vertices of an RC Tree.
  //If no edge (u,v) exists, throw an exception. 
  if(u > v) std::swap(u,v); 
  assert(edge_exists(u,v) && u != v && u >= 0 && v < n);
  RCCluster<aug_t>* edge_to_delete = nullptr;
  for(int i = 0; i < DEGREE_BOUND; i++){
    // Loop through adjacency lists of vertices u and v in the tree 
    // and delete the edge(u,v) between them.
    auto neighbor_cluster_u = contraction_tree[u][0][i];
    auto neighbor_cluster_v = contraction_tree[v][0][i];
    if(neighbor_cluster_u != nullptr && neighbor_cluster_u->bv_size() == 2){
      if((GET_NEIGHBOR(u, (*neighbor_cluster_u))) == v){
        edge_to_delete = neighbor_cluster_u;
        contraction_tree[u][0][i] = nullptr;
      }
    }

    if(neighbor_cluster_v != nullptr && neighbor_cluster_v->bv_size() == 2){
      auto val = GET_NEIGHBOR(v, (*neighbor_cluster_v));
      if(val == u){
        contraction_tree[v][0][i] = nullptr;
      }
    }
  }
  vec_insert(&affected, u);
  vec_insert(&affected, v); 

  vec_insert(&to_delete, edge_to_delete);
  roots.clear();
  update();
}

template<typename aug_t>
void RCTree<aug_t>::rake(vertex_t vertex, int round){
  /* An implementation of rake. Will take a degree 1 affected vertex and contract the vertex with its
 * adjacent neighbor by constructing a new cluster made of all unary clusters adjacent to the 
 * vertex and the single binary cluster between it and the neighbor it contracts into.*/

  //Get boundary vertex
  auto neighbors = contraction_tree[vertex][round];
  vertex_t neighbor; //To be initialized later, helpful for doing contractions.
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(id); //New augmented cluster to be used for replacement. 

  // Go through all neighboring clusters and aggregate values onto self.
  for(int i = 0; i < DEGREE_BOUND; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = new_cluster;
    if(is_edge(neighbors[i])){
      neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
    }
  }

  new_cluster->add_boundary(neighbor);
  new_cluster->rep_vertex = vertex;
  //Add to the neighbor this new unary cluster.

  add_neighbor(round + 1, new_cluster, vertex);
  delete representative_clusters[vertex];
  representative_clusters[vertex] = new_cluster; //need to free the old tuple

  // Add other vertex I am raking onto into affected set if 
  // the cluster of the vertex being raked in the next round
  // is not a unary cluster. Even if this vertex raked onto the exact
  // same vertex in the previous contraction the other vertex still needs
  // to be added to the affected set as the new cluster could have a different 
  // value.
  vec_insert(&affected, neighbor); 

  round_contracted[vertex] = round;
  // Free up extra memory used by sequence of vertex that contracted.
  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete[] contraction_tree[vertex][i];
  }
  contraction_tree[vertex].resize(round + 1);
  vec_erase(&affected, vertex);
}

template<typename aug_t>
void RCTree<aug_t>::compress(vertex_t vertex, int round){
  auto neighbors = contraction_tree[vertex][round];
  RCCluster<aug_t>* new_cluster = new RCCluster<aug_t>(id); //New augmented cluster to be used for replacement. 

  // Exactly the same as rake, but instead of a unary cluster being added to the list of 1 boundary vertex, it is a now 
  // a new binary cluster being added to the adjacency lists of both of its boundary_vertexes.

  // Go through all neighboring clusters and aggregate values onto self.
  for(int i = 0; i < DEGREE_BOUND; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = new_cluster;
    if(neighboring_cluster.bv_size() == 2){
      auto neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
      new_cluster->add_boundary(neighbor);
    }
  }

  new_cluster->rep_vertex = vertex;
  auto neighbor1 = new_cluster->boundary_vertexes[0];
  auto neighbor2 = new_cluster->boundary_vertexes[1];
  //Add to the neighbors this new binary cluster.
  add_neighbor(round + 1, new_cluster, vertex);
  delete representative_clusters[vertex];
  representative_clusters[vertex] = new_cluster;  
  
  vec_insert(&affected, neighbor1);
  vec_insert(&affected, neighbor2);

  round_contracted[vertex] = round;

  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete[] contraction_tree[vertex][i];
  }

  contraction_tree[vertex].resize(round + 1);
  vec_erase(&affected, vertex);
}

template<typename aug_t>
void RCTree<aug_t>::finalize(vertex_t vertex, int round){
  // Called when a vertex is the only uncontracted member left in the tree. Creates a new cluster for the
  // vertex which serves as the root of a particular RCTree in the forest.

  auto neighbors = contraction_tree[vertex][round];
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(id);  
  for(int i = 0; i < DEGREE_BOUND; i++){
    if(neighbors[i] == nullptr) continue;
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = new_cluster;
  }
  new_cluster->rep_vertex = vertex;
  roots.push_back(vertex);
  delete representative_clusters[vertex];
  representative_clusters[vertex] = new_cluster;
  round_contracted[vertex] = round;
  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete[] contraction_tree[vertex][i];
  }
  contraction_tree[vertex].resize(round + 1);
  vec_erase(&affected, vertex);
}

template<typename aug_t>
void RCTree<aug_t>::update() {
  // Based on set of affected vertices, create an MIS of the affected vertices, which induce
  // a subtree on the original tree and recontract them accordingly. Then determine new
  // set of affected vertices from vertices that contracted and did not contract and recurse.

  // Spread affection and redo contractions till there are no affected
  // vertices.

  int round = 0;

  while(!affected.empty()){

    //is_valid_induced_tree(round); - For testing purposes only
    spread_affection(round);
    MIS(round);
    //is_valid_MIS(round); - For testing purposes only
    
    for(auto vertex: affected){vec_insert(&to_delete, representative_clusters[vertex]);}

    #ifdef COLLECT_ROOT_CLUSTER_STATS
        if (rc_root_clusters_histogram.find(rc_root_clusters_histogram[round].size()) == rc_root_clusters_histogram.end())
            rc_root_clusters_histogram[rc_root_clusters_histogram[round].size()] = 1;
        else
            rc_root_clusters_histogram[rc_root_clusters_histogram[round].size()] += 1;
    #endif
    //Insert new neighbor list for next round for each affected vertex if it does not already
    //have a list inserted into the sequence to allow for insertions in its next round, as these are presently
    //the only ones that need to be updated
    for(auto vertex : affected){
      if(contraction_tree[vertex].size() < round + 2){
        contraction_tree[vertex].push_back(new RCCluster<aug_t>*[DEGREE_BOUND]{nullptr});
      }

      for(int i = 0; i < DEGREE_BOUND; i++){
        contraction_tree[vertex][round + 1][i] = nullptr;
        auto cluster = contraction_tree[vertex][round][i];
        if(is_edge(cluster)){
          auto neighbor = GET_NEIGHBOR(vertex, (*cluster));
          if(!vector_contains(&affected, neighbor) && !contracts(neighbor, round, 0)){
            clear_deleted_clusters(neighbor, round + 1);
          }
        }
      }
    }


    for(vertex_t v : affected) uncontracting.push_back(v); // Keeps track of vertices that did not contract in this round
    for(auto vertex : maximal_set){
      //Redo contractions for all vertices in maximal set.
      //Remove vertices from set of affected vertices.

      vec_erase(&uncontracting, vertex); // if vertex does not contract, neighbors need to update cluster associated with it. 
      auto vertex_degree = get_degree(vertex, round);
      if(vertex_degree == 1) {
        rake(vertex, round);
      } else if(vertex_degree == 2){ 
        compress(vertex, round);
      }else if(vertex_degree == 0) {
        // if the degree of a vertex is not 1 or 2 then it must be 0.
        // We have found the new root of the RC Tree and as such, we can
        // now finalize the vertex as the root of the new tree.
        finalize(vertex, round);
        break;
      }
    }

    // For each vertex that did not contract, add the binary cluster between
    // it and its neighbors to the adjacency lists of its neighbors in the
    // next round.
    for(auto vertex : uncontracting){
      for(int i = 0; i < DEGREE_BOUND; i++){
        auto cluster = contraction_tree[vertex][round][i];

        if(cluster == nullptr) continue; 

        if(cluster->bv_size() == 1) {
          add_neighbor(round + 1, cluster, vertex);
        }

        if(cluster->bv_size() == 2) {
          auto dereferenced_cluster = *cluster;
          auto neighbor = GET_NEIGHBOR(vertex, (*cluster));
          // If unaffected neighbor
          if(!vector_contains(&affected, neighbor) && !vector_contains(&maximal_set, neighbor)){
            if(contracts(neighbor, round ,1)){
              for(int i = 0; i < DEGREE_BOUND; i++){
                if(contraction_tree[vertex][round + 1][i] == nullptr){
                  contraction_tree[vertex][round + 1][i] = representative_clusters[neighbor];
                  break;
                }
              }
            } else {
              add_neighbor(round + 1, cluster, vertex);
              if(contracts(vertex, round, 0)){
                vec_insert(&affected, neighbor); 
              }
            }
          } else if (vector_contains(&affected, neighbor)){
            bool already_added = false;
            for(int i = 0; i < DEGREE_BOUND; i++){
              if(contraction_tree[neighbor][round + 1][i] == cluster){
                already_added = true;
                break;
              }
            }
            if(!already_added){
              add_neighbor(round + 1, cluster, vertex); 
            }
          }
        }
      }
    }
    round += 1;
    maximal_set.clear();
    uncontracting.clear();
  }
  to_delete.clear();
}

template<typename aug_t>
aug_t RCTree<aug_t>::path_query(vertex_t u, vertex_t v){
  assert(u != v);
  RCCluster<aug_t>* curr = round_contracted[u] < round_contracted[v] ? representative_clusters[v] : representative_clusters[u];
  RCCluster<aug_t>* curr2 = curr == representative_clusters[u]? representative_clusters[v] : representative_clusters[u];
  std::vector<RCCluster<aug_t>*> parents1;
  std::vector<RCCluster<aug_t>*> parents2;
  // Find LCA
  while(curr != nullptr){
    parents1.push_back(curr);
    curr = curr->parent;
  }
  while(curr2 != nullptr){ 
    parents2.push_back(curr2);
    curr2 = curr2->parent;
  } 
  int i = parents1.size() - 1, j = parents2.size() - 1;
  /*for(int i = 0; i < parents1.size(); i++){
    std::cout << parents1[i] << " ";
  }
  std::cout << "\n";
  for(int i = 0; i < parents2.size(); i++){
    std::cout << parents2[i] << " ";
  }
  std::cout << "\n";*/
  while(j > 0 && i > 0 && parents1[i-1] == parents2[j-1]){
    i-=1;
    j-=1;
  } 
  RCCluster<aug_t>* lca;
  if(i == 0){lca = parents1[i];} else {lca = parents2[j];}
  
  aug_t path_u1,path_u2, path_v1, path_v2; 
  path_u1 = path_u2 = path_v1 = path_v2 = id;
  // Compute initial cluster paths - base case of induction.
  for(int i = 0; i < DEGREE_BOUND;i++){
    auto cluster = contraction_tree[u][round_contracted[u]][i];
    if(cluster && cluster->bv_size() == 2){
      if(path_u1 == id){
        path_u1 = cluster->aug_val; 
      } else{
        path_u2 = cluster->aug_val;
      }
    }
  }

  for(int i = 0; i < DEGREE_BOUND;i++){
    auto cluster = contraction_tree[v][round_contracted[v]][i];
    if(cluster && cluster->bv_size() == 2){
      if(path_v1 == id){
        path_v1 = cluster->aug_val; 
      } else{
        path_v2 = cluster->aug_val;
      }
    }
  }
  
  curr = representative_clusters[u]; curr2 = representative_clusters[v];
  vertex_t boundary_u_1, boundary_u_2, boundary_v_1, boundary_v_2; 
  if(get_degree(u, round_contracted[u]) >= 1){
    boundary_u_1 = curr->boundary_vertexes[0];
    if(get_degree(u, round_contracted[u]) == 2) boundary_u_2 = curr->boundary_vertexes[1];
  }
  if(get_degree(v, round_contracted[v]) >= 1){
    boundary_v_1 = curr2->boundary_vertexes[0];
    if(get_degree(v, round_contracted[v]) == 2) boundary_v_2 = curr2->boundary_vertexes[1];
  } 
  
  
  //std::cout << lca << "\n";
  // Compute Path query.
  while((curr->parent != lca && curr != lca) || (curr2->parent != lca && curr2 != lca)){
    //std::cout << curr << " " << curr2 << "\n";
    if(curr->parent != lca && curr != lca){
      auto parent = curr->parent;
      vertex_t parent_id = parent->rep_vertex, curr_id = curr->rep_vertex;
      int degree_u = get_degree(curr_id, round_contracted[curr_id]), 
          degree_parent = get_degree(parent_id, round_contracted[parent_id]);
      
      // NOTE(ATHARVA): Make this all into one function.
      // Binary to Binary
      if(degree_u == 2 && degree_parent == 2){
        vertex_t incorrect_boundary = boundary_u_1 != parent_id ? boundary_u_1 : boundary_u_2;
        for(int i = 0; i < DEGREE_BOUND; i++) {
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            vertex_t neighbor = GET_NEIGHBOR(parent_id, (*cluster));
            if(neighbor != incorrect_boundary){
              if(incorrect_boundary != boundary_u_1){
                path_u1 = f(path_u1, cluster->aug_val);
                boundary_u_1 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
              } else{
                path_u2 = f(path_u2, cluster->aug_val);
                boundary_u_2 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
              }
            }
          }
        } 
        //Binary to Unary
      } else if(degree_u == 2 && degree_parent == 1){
        if(boundary_u_1 == parent_id){ 
          path_u1 = path_u2;
          boundary_u_1 = boundary_u_2;
        }
        path_u2 = id; boundary_u_2 = -1;
      } 
      // Unary to Binary
      else if(degree_u == 1 && degree_parent == 2){
        boundary_u_1 = parent->boundary_vertexes[0];
        boundary_u_2 = parent->boundary_vertexes[1];

        for(int i = 0; i < DEGREE_BOUND; i++){
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            vertex_t neighbor = GET_NEIGHBOR(parent_id, (*cluster));
            if(neighbor == boundary_u_1){
              path_u1 = f(path_u1, cluster->aug_val);
            } else {
              path_u2 = f(path_u2, cluster->aug_val);
            }
          }
        }
      } // Unary to Unary
      else {
        for(int i = 0; i < DEGREE_BOUND; i++){
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            path_u1 = f(path_u1, cluster->aug_val);
            boundary_u_1 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
          }
        } 
      }
      curr = parent;
    }

    if(curr2->parent != lca && curr2 != lca){ 
      auto parent = curr2->parent;
      vertex_t parent_id = parent->rep_vertex, curr_id = curr2->rep_vertex;
      int degree_v = get_degree(curr_id, round_contracted[curr_id]), degree_parent = get_degree(parent_id, round_contracted[parent_id]);

      // Binary to Binary
      if(degree_v == 2 && degree_parent == 2){
        vertex_t incorrect_boundary = boundary_v_1 != parent_id ? boundary_v_1 : boundary_v_2;
        for(int i = 0; i < DEGREE_BOUND; i++) {
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            vertex_t neighbor = GET_NEIGHBOR(parent_id, (*cluster));
            if(neighbor != incorrect_boundary){
              if(incorrect_boundary != boundary_v_1){
                path_v1 = f(path_v1, cluster->aug_val);
                boundary_v_1 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
              } else{
                path_v2 = f(path_v2, cluster->aug_val);
                boundary_v_2 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
              }
            }
          }
        } 
        //Binary to Unary
      } else if(degree_v == 2 && degree_parent == 1){ 
        if(boundary_v_1 == parent_id){ 
          path_v1 = path_v2;
          boundary_v_1 = boundary_v_2;
        }
        path_v2 = id; boundary_v_2 = -1;

      } 
      // Unary to Binary
      else if(degree_v == 1 && degree_parent == 2){
        boundary_v_1 = parent->boundary_vertexes[0];
        boundary_v_2 = parent->boundary_vertexes[1];

        for(int i = 0; i < DEGREE_BOUND; i++){
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            vertex_t neighbor = GET_NEIGHBOR(parent_id, (*cluster));
            if(neighbor == boundary_v_1){
              path_v1 = f(path_v1, cluster->aug_val);
            } else {
              path_v2 = f(path_v2, cluster->aug_val);
            }
          }
        }
      } // Unary to Unary
      else {
        for(int i = 0; i < DEGREE_BOUND; i++){
          auto cluster = contraction_tree[parent_id][round_contracted[parent_id]][i];
          if(cluster && cluster->bv_size() == 2){
            path_v1 = f(path_v1, cluster->aug_val);
            boundary_v_1 = cluster->boundary_vertexes[0] != parent_id ? cluster->boundary_vertexes[0] : cluster->boundary_vertexes[1];
          }
        } 
      }
      curr2 = parent;
    }
    if(curr->parent != lca || curr2->parent != lca){
      continue;
    }
  }
  aug_t final_path = id;
  if(boundary_u_1 == lca->rep_vertex){final_path = f(final_path, path_u1);} else if(boundary_u_2 == lca->rep_vertex){final_path = f(final_path, path_u2);}
  if(boundary_v_1 == lca->rep_vertex){final_path = f(final_path, path_v1);} else if(boundary_v_2 == lca->rep_vertex){final_path = f(final_path, path_v2);}
  return final_path;
};

template<typename aug_t>
int RCCluster<aug_t>::bv_size(){
  int count = 0;
  for(int i = 0; i < 2;++i){
    if(boundary_vertexes[i] != MAX_VERTEX_T) count++;
  }
  return count;
};

template<typename aug_t>
void RCCluster<aug_t>::add_boundary(vertex_t v){
  assert(boundary_vertexes[1] == MAX_VERTEX_T);
  if(boundary_vertexes[0] == MAX_VERTEX_T){boundary_vertexes[0] = v;} else{boundary_vertexes[1] = v;}
};

}
