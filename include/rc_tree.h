#include <cstdio>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <parlay/sequence.h>
#include <stdexcept>
#include <tuple>
#include <unordered_set>
#include<vector>
#include<parlay/monoid.h>

//TODO: other methods that deal with degree, neighbors, etc. need to be updated to deal with 
//PROCESSING flag and vertices that are no longer alive.

typedef int32_t vertex_t;

#define GET_NEIGHBOR(source, cluster) cluster.boundary_vertexes[0] != source ? cluster.boundary_vertexes[0] : cluster.boundary_vertexes[1]
#define PROCESSING -1 //Flag that indicates and affected vertex is still being processed and has not contracted yet.

template<typename aug_t>
struct RCCluster {
public:
  aug_t aug_val; 
  std::vector<int> boundary_vertexes;

  vertex_t parent;


  RCCluster(aug_t _aug_val){ aug_val = _aug_val; parent = -1;}
  RCCluster(aug_t _aug_val, vertex_t _parent){
    aug_val = _aug_val;
    _parent = parent;
  }
};

template<typename aug_t>
class RCTree {
public: 
  int degree_bound, n; //n is number of vertices in the tree.

  std::unordered_set<int>* affected = new std::unordered_set<int>();
  parlay::sequence<parlay::sequence<RCCluster<aug_t>**>> contraction_tree; // RCCluster pointer not int - Made public for testing, will actually be private.  
  std::tuple<RCCluster<aug_t>*,int>* representative_clusters; //Stores the representative cluster and round contracted for each vertex.

  void is_valid_MIS(std::unordered_set<int> maximal_set, int round);
  void add_neighbor(int round, RCCluster<aug_t>* cluster, vertex_t v);
  int get_degree(int v, int round);
  bool contracts(vertex_t v, int round, bool version);
  void spread_affection(int round);
  bool affected_by_dependence(int vertex, int round, std::unordered_set<int> *current_affected);
  std::unordered_set<int>* MIS(int round);
  int find_update_neighbor(vertex_t neighbor, vertex_t v, int round);
  void rake(vertex_t vertex, int round);
  void compress(vertex_t vertex, int round);
  void update();
  int neighbor_count(RCCluster<aug_t>** neighbors);
  void finalize(vertex_t v, int round);
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
  
  // Initialize round 0 adjacency list for each vertex in the tree.
  for(int i = 0; i < _n; i++){
    contraction_tree.push_back(parlay::sequence<RCCluster<aug_t>**>());
  }
  representative_clusters = (std::tuple<RCCluster<aug_t>*, int>*) calloc(_n, sizeof(std::tuple<RCCluster<aug_t>*, int>)); 
  for(int i = 0; i < n; i++){ 
    contraction_tree[i].push_back(new RCCluster<aug_t>*[3]{nullptr});
    representative_clusters[i] = std::tuple<RCCluster<aug_t>*, int>{(RCCluster<aug_t>*) new RCCluster<aug_t>(0), 0};
  }
}

/* ------- HELPER METHODS ------- */
template<typename aug_t>
int RCTree<aug_t>::get_degree(int v, int round) {
  // Get the degree of a vertex at a particular round i.e. no. of binary clusters the
  // vertex's adjacency list at a round.
  int degree = 0;
  for(int i = 0; i < degree_bound; i++){
    auto cluster = contraction_tree[v][round][i];
    if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) degree++;
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
    return std::get<1>(representative_clusters[v]) == round; 
  }
  return std::get<1>(representative_clusters[v]) <= round &&
  std::get<1>(representative_clusters[v]) != PROCESSING;
}

template<typename aug_t>
void RCTree<aug_t>::add_neighbor(int round, RCCluster<aug_t> *cluster, vertex_t contracted_v){
  //For a given RCCluster, add cluster to adjacency lists of boundary vertices of the cluster
  //in that round. This was written to be a one size fits all method and can essentially
  //allow the caller of the method to write a cluster at any level of the contraction tree.

  // At level 0, the method will try to add the cluster if and only if there is an empty
  // spot in the adjacency lists of its boundary vertices where the cluster can be added.
  // At any level above that it will work as follows:
  
  //The method will first attempt to find a cluster to replace in the adjacency list and will replace
  //any one of the following. For any adjacency list only one of the following could happen:
  //Clusters viable to be replaced are defined as follows: 
  // 1) If the method is trying to add the contracted cluster of a vertex, it will replace
  // the contracted cluster of that same vertex in the adjacency list.
  //
  // 2) If a cluster in the adjacency list has the same boundary vertices as the cluster we are trying to insert.
  //
  // 3) If either of the boundary vertices of a cluster in the adjacency list are the equivalent to the contracted vertex
  // the cluster we are attempting to add represents. - Necessary to allow unary clusters to replace binary clusters.
  //
  // 4) The last case is quite tricky. The idea is that a situation could arise in which the cluster of a contracting
  // vertex is attempting to replace another previously contracted representative cluster of one of its boundary vertices.
  // This occurs when for example in the case given by line 130 in rc_tree_tests.cpp.
  //
  // Lastly, if we don't find a cluster to replace, then there MUST exist an empty spot in the adjacency list 
  // to which the cluster is added.
  //
  // In any other case, this is an invalid call to the method and it throws an exception to indicate this.

  for(vertex_t boundary_v : cluster->boundary_vertexes){
    if(boundary_v == contracted_v) continue;
    auto adjacencies = contraction_tree[boundary_v][round]; //Get adjacency list of a boundary vertex of the cluster at a round.

    bool vertex_added = false; // Flag that keeps track of whether some vertex was replaced in adjacencies.
    int null_idx = -1; // To add to null index iff it couldn't be added to any other slots in the adjacencies list.
    for(int i = 0; i < degree_bound; i++){
      auto neighbor_cluster = adjacencies[i];
      if(neighbor_cluster == nullptr){
        null_idx = i;
      } else {
        // Needs a fix to account for unary clusters that need to overwrite binary clusters.

        auto dereferenced_cluster = *cluster; // Needs to be done to allow a call to GET_NEIGHBOR on the cluster.
        // Conditional to check for the cases for replacement indicated above.
        if((contracted_v != -1 &&
          std::get<0>(representative_clusters[contracted_v]) == neighbor_cluster)
          || (neighbor_cluster->boundary_vertexes.size() == 2
          && boundaries_equal(neighbor_cluster, cluster))
          || (neighbor_cluster->boundary_vertexes[0] == contracted_v ||
          neighbor_cluster->boundary_vertexes[1] == contracted_v)
          || (cluster->boundary_vertexes.size() == 2 && 
              neighbor_cluster == std::get<0>(representative_clusters[GET_NEIGHBOR(boundary_v, dereferenced_cluster)]))){
          adjacencies[i] = cluster;
          vertex_added = true; // Add cluster, change flag and break.
          break;
        }
      
      } 
    }
    if(!vertex_added && null_idx == -1){
      std::cout << contracted_v << "\n";
      throw std::invalid_argument("Vertex could not be added to the tree");
    }
    if(!vertex_added){
      adjacencies[null_idx] = cluster;
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
void RCTree<aug_t>::is_valid_MIS(std::unordered_set<int> maximal_set, int round){
  //Given a computed MIS, check to see if the MIS is maximal and no
  //vertices added are neighbors of each other.

  //Checks to see if the set is maximal by going through all vertexes in the
  //affected set and checking to see if atleast one neighbor of every vertex is in
  //the maximal independent set.
  for(auto vertex : *affected){
    if(contraction_tree[vertex][round] != nullptr){ 
      bool neighbor_in_max = false;
      for(int i = 0; i < degree_bound; i++){
        if(contraction_tree[vertex][round][i] != nullptr && 
          contraction_tree[vertex][round][i]->boundary_vertexes.size() == 2){
          auto dereferenced_cluster = *contraction_tree[vertex][round][i];
          int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
          if(maximal_set.count(vertex) == 0 && 
            ((contracts(neighbor, round, 0)  && affected->count(neighbor) == 0)|| 
            maximal_set.count(neighbor) == 1)){
            neighbor_in_max = true; // If either the vertex is next to an contracting unaffected neighbor or 
            // neighbor included in the MIS then increment count as this vertex is valid.
          }

          // If 2 adjacent neighbors have been added then also throw an exception.
          // as the set is not independent.
          if(maximal_set.count(vertex) == 1 && 
            maximal_set.count(neighbor) != 0){
            throw std::invalid_argument("2 adjacent vertices have been added to the MIS");
          }
        }
      }

      //If no neighbors are in the maximal set then this vertex should have been
      //in the maximal set. Thus the MIS is invalid
      if(maximal_set.count(vertex) == 0 && get_degree(vertex, round) < 3 &&
          !neighbor_in_max){
        throw std::invalid_argument("The set is not maximal :-(");
      }
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
  std::unordered_set<int> *new_affected = new std::unordered_set<int>(*affected);

  for(auto vertex : *affected) {
    //Iterate through all adjacencies of affected vertex
    for(int i = 0; i < degree_bound; i++) {

      RCCluster<aug_t>* neighboring_cluster = contraction_tree[vertex][round][i];

      if(neighboring_cluster != nullptr && neighboring_cluster->boundary_vertexes.size() == 2) {
        RCCluster<aug_t> dereferenced_cluster = *neighboring_cluster; // Why can't I just pass in *cluster directly to GET_NEIGHBOR?
        int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster); // Neighbor of the affected vertex.
        //Check to see if neighbor is already affected and if all contracting neighbors are affected.
        if(new_affected->count(neighbor) == 0 && affected_by_dependence(neighbor, round, new_affected)){ 
          new_affected->insert(neighbor);
          representative_clusters[neighbor] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[neighbor]), PROCESSING};
          contraction_tree[neighbor].resize(round + 1);
        }
      } 
    }
  }
  delete affected;
  affected = new_affected;
}

template<typename aug_t>
bool RCTree<aug_t>::affected_by_dependence(int vertex, int round, std::unordered_set<int> *current_affected){
  /* Determine if a vertex becomes affected as result of contracting neighbors added to the affected set.
   * Return true if the vertex is now affected, false otherwise.*/ 
  if(contracts(vertex, round, 0)){
    return false;
  }
  bool spreads = true;
  for(int i = 0; i < degree_bound; i++){ //Neighbors of the neighbor of the affected vertex.
    auto neighbor_cluster = contraction_tree[vertex][round][i]; 
    if(neighbor_cluster != nullptr){
      RCCluster<aug_t> dereferenced_cluster = *neighbor_cluster;
      int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster); // Neighbor of the affected vertex.
      // If my neighbor does not contract
      if(std::get<0>(representative_clusters[neighbor])->parent == vertex && !current_affected->count(neighbor)) {
        spreads = false;
        break;
      }
    }
  }
  return spreads;
}

template<typename aug_t>
std::unordered_set<int>* RCTree<aug_t>::MIS(int round){
  //pick a random vertex
  //Include vertex in MIS
  //Delete Neighbors
  //Return

  std::unordered_set<int>* maximal_set = new std::unordered_set<int>(*affected);
  for(auto random_vertex : *affected){
    // If the random vertex is not in the maximal_set then the vertex has already
    // been removed as a neighbor of some other vertex and we can't do any ops on it.
    if(maximal_set->count(random_vertex) == 0) continue;
    // A few cases to deal with :
    // 1) If the vertex is next to a contracting unaffected vertex it is not to be considered.
    // 2) If the vertex cannot contract, it must remain in the set of affected vertices.

    // Case if next to unaffected contracting vertex.
    bool adjacent_unaffected_contracting = false;
    for(int i = 0; i < degree_bound; i++){

      auto cluster = contraction_tree[random_vertex][round][i];
      if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) {
        auto dereferenced_cluster = *cluster;
        int neighbor = GET_NEIGHBOR(random_vertex, dereferenced_cluster); 

        //Neighbor not in MIS and neighbor contracts in this round.
        if(!maximal_set->count(neighbor) && contracts(neighbor, round, 0)){
          maximal_set->erase(random_vertex); // Remove from MIS, if so.
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
        if(cluster != nullptr && cluster->boundary_vertexes.size() == 2){
          auto dereferenced_cluster = *cluster;
          int neighbor = GET_NEIGHBOR(random_vertex, dereferenced_cluster);
          if(maximal_set->count(neighbor)) maximal_set->erase(neighbor);
        }
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

  // Create RC cluster to represent edge at level 0 between u and v.
  RCCluster<aug_t>* new_edge = (RCCluster<aug_t>*) new RCCluster<aug_t>(weight);
  new_edge->boundary_vertexes.push_back(u);
  new_edge->boundary_vertexes.push_back(v);

  add_neighbor(0, new_edge, -1);

  affected->insert(u);          // Insert initial affected vertices u and v.
  affected->insert(v);

  
  // Change representative_clusters tuples of each vertex to indicate the vertex is being processed and as
  // such the round in which it contracts now has not yet been computed.
  representative_clusters[u] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[u]), PROCESSING};
  representative_clusters[v] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[v]), PROCESSING};

  update(); 
}

template<typename aug_t>
void RCTree<aug_t>::cut(int u, int v) {
  //Delete an edge between 2 vertices of an RC Tree.
  //If no edge (u,v) exists, throw an exception.
  for(int i = 0; i < degree_bound; i++){
    auto neighbor_cluster_u = contraction_tree[u][0][i];
    auto neighbor_cluster_v = contraction_tree[v][0][i];
    if(neighbor_cluster_u != nullptr && neighbor_cluster_u->boundary_vertexes.size() == 2){
      auto dereferenced_cluster = *neighbor_cluster_u; 
      if(GET_NEIGHBOR(u, dereferenced_cluster) == v){
        contraction_tree[u][0][i] = nullptr;
      }
    }

    if(neighbor_cluster_v != nullptr && neighbor_cluster_v->boundary_vertexes.size() == 2){
      auto dereferenced_cluster = *neighbor_cluster_v; 
      if(GET_NEIGHBOR(v, dereferenced_cluster) == u){
        contraction_tree[v][0][i] = nullptr;
      }
    }
  }
  affected->insert(u);
  affected->insert(v);

  representative_clusters[u] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[u]), PROCESSING};
  representative_clusters[v] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[v]), PROCESSING};

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
    new_cluster->aug_val += neighbors[i]->aug_val; // Needs to be a associative binary op passed in by the user.
    neighbors[i]->parent = vertex;
    if(neighboring_cluster.boundary_vertexes.size() == 2){
      neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
    }
  }

  new_cluster->boundary_vertexes.push_back(neighbor);
  //Add to the neighbor this new unary cluster.

  add_neighbor(round + 1, new_cluster, vertex);
  representative_clusters[vertex] = std::tuple{new_cluster, round}; //need to free the old tuple

  // Add other vertex I am raking onto into affected set if 
  // the cluster of the vertex being raked in the next round
  // is not a unary cluster. Even if this vertex raked onto the exact
  // same vertex in the previous contraction the other vertex still needs
  // to be added to the affected set as the new cluster could have a different 
  // value.

  affected->insert(neighbor);

  // Free up extra memory used by sequence of vertex that contracted.
  for(int i = round + 1; i < contraction_tree[vertex].size(); i++){
    delete contraction_tree[vertex][i];
  }
  contraction_tree[vertex].resize(round + 1);
  representative_clusters[neighbor] = std::tuple{std::get<0>(representative_clusters[neighbor]), PROCESSING};
  affected->erase(vertex);
}

template<typename aug_t>
void RCTree<aug_t>::compress(vertex_t vertex, int round){
  auto neighbors = contraction_tree[vertex][round];
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(0); //New augmented cluster to be used for replacement. - Monoid struct in parlay include.

  // Exactly the same as rake, but instead of a unary cluster being added to the list of 1 boundary vertex, it is a now 
  // a new binary cluster being added to the adjacency lists of both of its boundary_vertexes.

  // Go through all neighboring clusters and aggregate values onto self.
  for(int i = 0; i < degree_bound; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val += neighboring_cluster.aug_val; // Needs to be a associative binary op passed in by the user.
    neighboring_cluster.parent = vertex;
    if(neighboring_cluster.boundary_vertexes.size() == 2){
      auto neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
      new_cluster->boundary_vertexes.push_back(neighbor);
    }
  }

  auto neighbor1 = new_cluster->boundary_vertexes[0];
  auto neighbor2 = new_cluster->boundary_vertexes[1];
  //Add to the neighbors this new binary cluster.
  add_neighbor(round + 1, new_cluster, vertex);
  representative_clusters[vertex] = std::tuple{new_cluster, round}; //need to free the old tuple

  affected->insert(neighbor1);
  affected->insert(neighbor2);
  contraction_tree[vertex].resize(round + 1);
  representative_clusters[neighbor1] = std::tuple{std::get<0>(representative_clusters[neighbor1]), PROCESSING};
  representative_clusters[neighbor2] = std::tuple{std::get<0>(representative_clusters[neighbor2]), PROCESSING};
  affected->erase(vertex);
}

template<typename aug_t>
void RCTree<aug_t>::finalize(vertex_t v, int round){
  // Called when a vertex is the only uncontracted member left in the tree. Creates a new cluster for the
  // vertex and adds to it to the representative_clusters map.
  auto neighbors = contraction_tree[v][round];
  RCCluster<aug_t> *new_cluster = new RCCluster<aug_t>(0, -1);  
  for(int i = 0; i < degree_bound; i++){
    if(neighbors[i] == nullptr) continue;
    auto neighboring_cluster = *neighbors[i];
    new_cluster->aug_val += neighbors[i]->aug_val; // Needs to be a associative binary op passed in by the user.
    neighboring_cluster.parent = v;
  }

  representative_clusters[v] = *(new std::tuple<RCCluster<aug_t>*, int>{new_cluster, round});
  contraction_tree[v].resize(round + 1);
  affected->erase(v);
}
template<typename aug_t>
void RCTree<aug_t>::update() {
  // Based on set of affected vertices,.count an MIS of the affected vertices, which induce
  // a subtree on the original tree and recontract them accordingly. Then determine new
  // set of affected vertices and recurse.

  // Spread affection and redo contractions till there are affected
  // vertices.

  int round = 0;


  while(!affected->empty()){

    spread_affection(round);
    std::unordered_set<int> *maximal_set = MIS(round);
    is_valid_MIS(*maximal_set, round);
    //Insert new neighbor list for next round for each affected vertex if it does not already
    //have a list inserted into the sequence to allow for insertions in its next round, as these are presently
    //the only ones that need to be updated.
    for(auto vertex : *affected){
      if(contraction_tree[vertex].size() < round + 2){
        contraction_tree[vertex].push_back(new RCCluster<aug_t>*[degree_bound]{nullptr});
        //Copy over neighbor list from previous round for each vertex.
        for(int i = 0; i < degree_bound; i++){
          contraction_tree[vertex][round + 1][i] = contraction_tree[vertex][round][i];

          auto prev_cluster = contraction_tree[vertex][round][i];
          if(prev_cluster == nullptr || prev_cluster->boundary_vertexes.size() == 1){
            continue;
          }
          auto dereferenced_cluster = *prev_cluster;
          auto neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
          if(contracts(neighbor,round,1) && affected->count(neighbor) == 0){
            contraction_tree[vertex][round + 1][i] = std::get<0>(representative_clusters[neighbor]); 
          }
        }
      } else {
        // If the sequence for the affected vertex already as an adjacency list for the next round
        // this adjacency list could contain clusters that do not exist anymore as a result
        // of contractions that could have happened. Delete these clusters by checking to see
        // if the boundary vertex they are adjacent to can be found as either the representative
        // vertex or boundary vertex of another cluster in the previous round. If not found, set this
        // entry in the adjacency list for the next round to be null.
        for(int i = 0; i < degree_bound;i++){
          if(contraction_tree[vertex][round + 1][i] != nullptr){
            auto cluster = contraction_tree[vertex][round][i];
            bool elem_found = false;
            for(int j = 0; j < degree_bound; j++){
              if(contraction_tree[vertex][round][j] != nullptr){
                auto dereferenced_cluster = *contraction_tree[vertex][round][j];
                auto neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);

                if(cluster == contraction_tree[vertex][round][j] || 
                  std::get<0>(representative_clusters[neighbor]) == cluster){
                  elem_found = true;
                  break;
                }
              }
            }

            if(!elem_found){
              contraction_tree[vertex][round][i] = nullptr;
            }
          }
        }
      }
    }

    std::unordered_set<int> uncontracting(*affected); // Keeps track of vertices that did not contract in this round
    for(auto vertex : *maximal_set){
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
        if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) {
          auto dereferenced_cluster = *cluster;
          auto neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
          if(contracts(neighbor, round, 1)) continue;

          add_neighbor(round + 1, cluster, vertex);

          affected->insert(neighbor);
          representative_clusters[neighbor] = std::tuple{std::get<0>(representative_clusters[neighbor]), PROCESSING};
        }
      }
    }

    round += 1;
    delete maximal_set;
  }
}
