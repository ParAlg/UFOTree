#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <parlay/sequence.h>
#include <stdexcept>
#include <unordered_set>
#include<vector>
#include "types.h"
#include "util.h"

#define GET_NEIGHBOR(source, cluster) cluster.boundary_vertexes[0] != source ? cluster.boundary_vertexes[0] : cluster.boundary_vertexes[1]

#define DEGREE_BOUND 3 // CHANGE THIS IN CASE YOU WANT A DIFFERENT DEGREE BOUND FOR YOUR PURPOSES.

template<typename aug_t>
struct RCCluster {
public:
  aug_t aug_val; 
  std::vector<vertex_t> boundary_vertexes;

  vertex_t parent;
  vertex_t rep_vertex;

  RCCluster(aug_t _aug_val){ aug_val = _aug_val; parent = MAX_VERTEX_T; rep_vertex = MAX_VERTEX_T;}
  RCCluster(aug_t _aug_val, vertex_t _parent){
    aug_val = _aug_val;
    parent = _parent;
    rep_vertex = MAX_VERTEX_T;
  }
  
  bool operator==(const RCCluster<aug_t> &other){
    return (aug_val == other.aug_val && 
            boundaries_equal(this, &other) &&
            parent == other.parent && rep_vertex == other.rep_vertex);
  }
};

template<typename aug_t>
class RCTree {
public:
  // Data structure fields.
  int degree_bound;
  vertex_t n; //n is number of vertices in the tree.
  aug_t id, d_val; //identity element, default value
  std::function<aug_t(aug_t, aug_t)> f; // Associative function
  std::vector<vertex_t> roots; //for debugging purposes only.
  std::unordered_set<vertex_t> affected;
  parlay::sequence<parlay::sequence<RCCluster<aug_t>**>> contraction_tree; 
  std::vector<RCCluster<aug_t>*> representative_clusters;   
  std::unordered_set<RCCluster<aug_t>*> to_delete;
  std::unordered_set<vertex_t> maximal_set;
  std::vector<int> round_contracted; 
  size_t max_space; // For space benchmarks.

  // Functions
  void is_valid_MIS(int round);
  void add_neighbor(int round, RCCluster<aug_t>* cluster, vertex_t v);
  int get_degree(vertex_t v, int round);
  bool contracts(vertex_t v, int round, bool version);
  void spread_affection(int round);
  bool affected_by_dependence(vertex_t vertex, int round, std::unordered_set<vertex_t> *current_affected);
  void MIS(int round);
  void rake(vertex_t vertex, int round);
  void compress(vertex_t vertex, int round);
  void update();
  int neighbor_count(RCCluster<aug_t>** neighbors);
  void finalize(vertex_t vertex, int round);
  vertex_t get_root(RCCluster<aug_t>* cluster);
  void clear_deleted_clusters(vertex_t vertex, int round);
  void is_valid_induced_tree(int round);
  bool connected(vertex_t u, vertex_t v);
  bool is_edge(RCCluster<aug_t> *cluster);
  bool edge_exists(vertex_t u, vertex_t v);
  RCTree<aug_t>** get_neighbors(vertex_t v);
  size_t space();
  /*RCTree::RCTree(vector<int[3]> tree, int n, int degree_bound); */
  // These are only ones that should really be public.
  RCTree(int _n, QueryType q = PATH, 
         std::function<aug_t(aug_t, aug_t)> f = [] (aug_t x, aug_t y) {return x + y;}, 
         aug_t id = 0, aug_t d_val = 0);
  void link(vertex_t u, vertex_t v, aug_t weight = 1);
  void cut(vertex_t u, vertex_t v);
};

template<typename aug_t>
// Constructor. 
RCTree<aug_t>::RCTree(int _n, QueryType _q, 
                      std::function<aug_t(aug_t, aug_t)> _f,
                      aug_t _id, aug_t _d_val) { 

  degree_bound = DEGREE_BOUND; n = _n; f = _f; id = _id; d_val = _d_val; 
  // Initialize round 0 adjacency list for each vertex in the tree.
  for(int i = 0; i < _n; i++){
    contraction_tree.push_back(parlay::sequence<RCCluster<aug_t>**>());
  }
  representative_clusters.resize(n); 
  round_contracted.resize(n); 
  for(int i = 0; i < n; i++){ 
    contraction_tree[i].push_back(new RCCluster<aug_t>*[degree_bound]{nullptr});
    representative_clusters[i] = (RCCluster<aug_t>*) new RCCluster<aug_t>(d_val);
    representative_clusters[i]->rep_vertex = i;
    round_contracted[i] = 0;
  }
}

/* ------- HELPER METHODS ------- */
template<typename aug_t>
int RCTree<aug_t>::get_degree(vertex_t v, int round) {
  // Get the degree of a vertex at a particular round i.e. no. of binary clusters in the
  // vertex's adjacency list at a round.
  int degree = 0;
  for(int i = 0; i < degree_bound; i++){
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
void RCTree<aug_t>::add_neighbor(int round, RCCluster<aug_t> *cluster, vertex_t contracted_v){
  /* Add a cluster into the adjacency list of its boundary vertex(es) at a particular
   * round*/
  for(vertex_t boundary_v : cluster->boundary_vertexes){
    //Get adjacency list of a boundary vertex of the cluster at a round.
    auto adjacencies = contraction_tree[boundary_v][round]; 
    auto vertex_added = false;

    for(int i = 0; i < degree_bound; i++){
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
vertex_t RCTree<aug_t>::get_root(RCCluster<aug_t>* cluster){
  // Returns the root of a cluster in the RC Tree if it exists.
  auto curr = cluster;
  while(curr->parent != MAX_VERTEX_T){
    curr = representative_clusters[curr->parent];
  }
  return curr->rep_vertex;
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
  return cluster != nullptr && cluster->boundary_vertexes.size() == 2;
}

template<typename aug_t>
bool RCTree<aug_t>::edge_exists(vertex_t u, vertex_t v){
  auto u_adjacencies = contraction_tree[u][0];
  for(int i = 0; i < degree_bound; i++){
    if(is_edge(u_adjacencies[i])){
      if((GET_NEIGHBOR(u, (*u_adjacencies[i]))) == v){
        bool contains_cluster = false;
        for(int j = 0; j < degree_bound; j++){
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
  for(int i = 0; i < degree_bound; i++){
    if(contraction_tree[vertex][round][i] != nullptr && 
        to_delete.count(contraction_tree[vertex][round][i]) == 1){
      contraction_tree[vertex][round][i] = nullptr;
    }

    if(is_edge(contraction_tree[vertex][round][i])){
      auto neighbor = GET_NEIGHBOR(vertex, (*contraction_tree[vertex][round][i]));
      if(affected.count(neighbor) == 1){
        contraction_tree[vertex][round][i] = nullptr;
      }
    }
  }
}
template<typename aug_t>
RCTree<aug_t>** RCTree<aug_t>::get_neighbors(vertex_t v){
  return contraction_tree[v][0];
}
/* ------------------------------- */

template<typename aug_t>
void RCTree<aug_t>::spread_affection(int round) {
  /* Determine the vertices we spread affection to via dependence as a result of changes
   * to the input and add them to the affected set, prior to building an MIS of the
   * affected vertices.
   */

  //Vertices that affection is spread to from orginal affected vertices
  std::unordered_set<vertex_t> new_affected(affected);

  for(auto vertex : affected) {
    //Iterate through all adjacencies of affected vertex
    for(int i = 0; i < degree_bound; i++) {

      RCCluster<aug_t>* neighboring_cluster = contraction_tree[vertex][round][i];

      if(is_edge(neighboring_cluster)) {
        int neighbor = GET_NEIGHBOR(vertex, (*neighboring_cluster)); // Neighbor of the affected vertex.
        //Check to see if neighbor is already affected and if all contracting neighbors are affected.
        if(new_affected.count(neighbor) == 0 && affected_by_dependence(neighbor, round, &new_affected)){ 
          new_affected.insert(neighbor);
        }
      } 
    }
  }
  for(vertex_t v : new_affected) affected.insert(v);
}

template<typename aug_t>
bool RCTree<aug_t>::affected_by_dependence(vertex_t vertex, int round, std::unordered_set<vertex_t> *current_affected){
  /* Determine if a vertex becomes affected as result of contracting neighbors added to the affected set.
   * Return true if the vertex is now affected, false otherwise.*/ 
  if(contracts(vertex, round, 0)) return false;

  bool spreads = true;
  for(int i = 0; i < degree_bound; i++){ //Neighbors of the neighbor of the affected vertex.
    auto neighbor_cluster = contraction_tree[vertex][round][i]; 
    if(is_edge(neighbor_cluster)){
      int neighbor = GET_NEIGHBOR(vertex, (*neighbor_cluster)); // Neighbor of the affected vertex.
      // If my neighbor does not contract
      if(contracts(neighbor, round, 0) && current_affected->count(neighbor) == 0) {
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

  for(auto v : affected) maximal_set.insert(v);

  for(auto random_vertex : affected){
    // If the random vertex is not in the maximal_set then the vertex has already
    // been removed as a neighbor of some other vertex and we can't do any ops on it.
    if(maximal_set.count(random_vertex) == 0) continue;
    // A few cases to deal with :
    // 1) If the vertex is next to a contracting unaffected vertex it is not to be considered.
    // 2) If the vertex cannot contract, it must remain in the set of affected vertices.

    // Case if next to unaffected contracting vertex.
    bool adjacent_unaffected_contracting = false;
    for(int i = 0; i < degree_bound; i++){
      auto cluster = contraction_tree[random_vertex][round][i];
      if(is_edge(cluster)) {
        int neighbor = GET_NEIGHBOR(random_vertex, (*cluster)); 

        //Neighbor not in affected set and neighbor contracts in this round.
        if(contracts(neighbor, round, 0) && affected.count(neighbor) == 0){
          maximal_set.erase(random_vertex); // Remove from MIS, if so.
          adjacent_unaffected_contracting = true;
          break;
        }
      }
    }

    // If not case 1 and vertex degree < 3, then this vertex is part of indep. set. 
    // Remove all neighbors.
    if(!adjacent_unaffected_contracting && get_degree(random_vertex, round) != 3){
      for(int i = 0; i < degree_bound; i++){
        auto cluster = contraction_tree[random_vertex][round][i];
        if(is_edge(cluster)){
          int neighbor = GET_NEIGHBOR(random_vertex, (*cluster));
          if(maximal_set.count(neighbor)) maximal_set.erase(neighbor);
        }
      }
    } else {
      maximal_set.erase(random_vertex);
    }
  }
}

template<typename aug_t>
void RCTree<aug_t>::link(vertex_t u, vertex_t v, aug_t weight) {
  // Add an edge between 2 Trees of the RC forest.
  assert(u >= 0 && v < n && !connected(u,v));

  // Create RC cluster to represent edge at level 0 between u and v.
  RCCluster<aug_t>* new_edge = (RCCluster<aug_t>*) new RCCluster<aug_t>(weight);
  new_edge->boundary_vertexes.push_back(u);
  new_edge->boundary_vertexes.push_back(v);

  add_neighbor(0, new_edge, MAX_VERTEX_T);
  affected.insert(u);          // Insert initial affected vertices u and v.
  affected.insert(v); 
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
  for(int i = 0; i < degree_bound; i++){
    // Loop through adjacency lists of vertices u and v in the tree 
    // and delete the edge(u,v) between them.
    auto neighbor_cluster_u = contraction_tree[u][0][i];
    auto neighbor_cluster_v = contraction_tree[v][0][i];
    if(neighbor_cluster_u != nullptr && neighbor_cluster_u->boundary_vertexes.size() == 2){
      if((GET_NEIGHBOR(u, (*neighbor_cluster_u))) == v){
        edge_to_delete = neighbor_cluster_u;
        contraction_tree[u][0][i] = nullptr;
      }
    }

    if(neighbor_cluster_v != nullptr && neighbor_cluster_v->boundary_vertexes.size() == 2){
      auto val = GET_NEIGHBOR(v, (*neighbor_cluster_v));
      if(val == u){
        contraction_tree[v][0][i] = nullptr;
      }
    }
  }
  affected.insert(u);
  affected.insert(v);
  to_delete.insert(edge_to_delete);
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
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(0); //New augmented cluster to be used for replacement. - Monoid struct in parlay include.

  // Go through all neighboring clusters and aggregate values onto self.
  for(int i = 0; i < degree_bound; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = vertex;
    if(is_edge(neighbors[i])){
      neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
    }
  }

  new_cluster->boundary_vertexes.push_back(neighbor);
  new_cluster->rep_vertex = vertex;
  //Add to the neighbor this new unary cluster.

  add_neighbor(round + 1, new_cluster, vertex);
  representative_clusters[vertex] = new_cluster; //need to free the old tuple

  // Add other vertex I am raking onto into affected set if 
  // the cluster of the vertex being raked in the next round
  // is not a unary cluster. Even if this vertex raked onto the exact
  // same vertex in the previous contraction the other vertex still needs
  // to be added to the affected set as the new cluster could have a different 
  // value.
  affected.insert(neighbor);

  round_contracted[vertex] = round;
  // Free up extra memory used by sequence of vertex that contracted.
  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete contraction_tree[vertex][i];
  }
  contraction_tree[vertex].resize(round + 1);
  affected.erase(vertex);
}

template<typename aug_t>
void RCTree<aug_t>::compress(vertex_t vertex, int round){
  auto neighbors = contraction_tree[vertex][round];
  RCCluster<aug_t>* new_cluster = new RCCluster<aug_t>(0); //New augmented cluster to be used for replacement. 

  // Exactly the same as rake, but instead of a unary cluster being added to the list of 1 boundary vertex, it is a now 
  // a new binary cluster being added to the adjacency lists of both of its boundary_vertexes.

  // Go through all neighboring clusters and aggregate values onto self.
  for(int i = 0; i < degree_bound; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = vertex;
    if(neighboring_cluster.boundary_vertexes.size() == 2){
      auto neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
      new_cluster->boundary_vertexes.push_back(neighbor);
    }
  }

  new_cluster->rep_vertex = vertex;
  auto neighbor1 = new_cluster->boundary_vertexes[0];
  auto neighbor2 = new_cluster->boundary_vertexes[1];
  //Add to the neighbors this new binary cluster.
  add_neighbor(round + 1, new_cluster, vertex);
  representative_clusters[vertex] = new_cluster; 

  affected.insert(neighbor1);
  affected.insert(neighbor2);

  round_contracted[vertex] = round;

  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete contraction_tree[vertex][i];
  }

  contraction_tree[vertex].resize(round + 1);
  affected.erase(vertex);
}

template<typename aug_t>
void RCTree<aug_t>::finalize(vertex_t vertex, int round){
  // Called when a vertex is the only uncontracted member left in the tree. Creates a new cluster for the
  // vertex which serves as the root of a particular RCTree in the forest.

  auto neighbors = contraction_tree[vertex][round];
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(0, MAX_VERTEX_T);  
  for(int i = 0; i < degree_bound; i++){
    if(neighbors[i] == nullptr) continue;
    new_cluster->aug_val = f(new_cluster->aug_val, neighbors[i]->aug_val); // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = vertex;
  }
  new_cluster->rep_vertex = vertex;
  roots.push_back(vertex);
  representative_clusters[vertex] = new_cluster;
  round_contracted[vertex] = round;
  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete contraction_tree[vertex][i];
  }
  contraction_tree[vertex].resize(round + 1);
  affected.erase(vertex);
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

    //is_valid_induced_tree(round);
    spread_affection(round);
    MIS(round);
    //is_valid_MIS(round);
    
    for(auto vertex: affected){to_delete.insert(representative_clusters[vertex]);}

    //Insert new neighbor list for next round for each affected vertex if it does not already
    //have a list inserted into the sequence to allow for insertions in its next round, as these are presently
    //the only ones that need to be updated
    for(auto vertex : affected){
      if(contraction_tree[vertex].size() < round + 2){
        contraction_tree[vertex].push_back(new RCCluster<aug_t>*[degree_bound]{nullptr});
      }

      for(int i = 0; i < degree_bound; i++){
        contraction_tree[vertex][round + 1][i] = nullptr;
        auto cluster = contraction_tree[vertex][round][i];
        if(is_edge(cluster)){
          auto neighbor = GET_NEIGHBOR(vertex, (*cluster));
          if(affected.count(neighbor) == 0 && !contracts(neighbor, round, 0)){
            clear_deleted_clusters(neighbor, round + 1);
          }
        }
      }
    }


    std::unordered_set<vertex_t> uncontracting(affected); // Keeps track of vertices that did not contract in this round
    for(auto vertex : maximal_set){
      //Redo contractions for all vertices in maximal set.
      //Remove vertices from set of affected vertices.

      uncontracting.erase(vertex); // if vertex does not contract, neighbors need to update cluster associated with it. 
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
      for(int i = 0; i < degree_bound; i++){
        auto cluster = contraction_tree[vertex][round][i];

        if(cluster == nullptr) continue; 

        if(cluster->boundary_vertexes.size() == 1) {
          add_neighbor(round + 1, cluster, vertex);
        }

        if(cluster->boundary_vertexes.size() == 2) {
          auto dereferenced_cluster = *cluster;
          auto neighbor = GET_NEIGHBOR(vertex, (*cluster));
          // If unaffected neighbor
          if(!affected.count(neighbor) && maximal_set.count(neighbor) == 0){
            if(contracts(neighbor, round ,1)){
              for(int i = 0; i < degree_bound; i++){
                if(contraction_tree[vertex][round + 1][i] == nullptr){
                  contraction_tree[vertex][round + 1][i] = representative_clusters[neighbor];
                  break;
                }
              }
            } else {
              add_neighbor(round + 1, cluster, vertex);
              if(contracts(vertex, round, 0)){
                affected.insert(neighbor);
              }
            }
          } else if (affected.count(neighbor) == 1){
            bool already_added = false;
            for(int i = 0; i < degree_bound; i++){
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
  }
  to_delete.clear();
}
