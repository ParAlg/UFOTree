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

  parlay::sequence<RCCluster<aug_t>> clusters;
  parlay::sequence<RCCluster<aug_t>> leaf_clusters;
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
  // Get the degree of a vertex at a particular round. 

  int degree = 0;
  for(int i = 0; i < degree_bound; i++){
    auto cluster = contraction_tree[v][round][i];
    if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) degree++;
  }

  return degree;
}
/*
template<typename aug_t>
int RCTree<aug_t>::neighbor_count(RCCluster<aug_t>** neighbors){
  int count = 0;
  for(int i = 0; i < degree_bound; i++){
    auto neighbor_ptr = neighbors[i];
    if(neighbor_ptr != nullptr && 
       neighbor_ptr->round_contracted != PROCESSING && neighbor_ptr->round_contracted < round){
      count += 1;
    }
  }
  return count;
}
*/
template<typename aug_t>
bool RCTree<aug_t>::contracts(vertex_t v, int round, bool version) {
  // Take in vertex and round number to determine if the vertex and which variant
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
  //For a given RCCluster, add cluster to neighbor lists of boundary vertices of that 
  //vertex at a particular round if neighbors is full, throw an exception.

  // If the cluster is not represented yet by any vertex v, passing in -1 for the v param
  // will allow the cluster to be added to a nullptr location of the adjacencies.

  for(vertex_t boundary_v : cluster->boundary_vertexes){
    if(boundary_v == contracted_v) continue;
    auto adjacencies = contraction_tree[boundary_v][round]; //Get neighbor list of boundary at round.

    bool vertex_added = false;
    int null_idx = -1; // To add to null index iff it couldn't be added to any other slots in the adjacencies list.
    for(int i = 0; i < degree_bound; i++){
      auto neighbor_cluster = adjacencies[i];
      if(neighbor_cluster == nullptr){
        null_idx = i;
      } else {
        // Needs a fix to account for unary clusters that need to overwrite binary clusters.

        
        if((contracted_v != -1 &&
          std::get<0>(representative_clusters[contracted_v]) == neighbor_cluster)
          || (neighbor_cluster->boundary_vertexes.size() == 2
          && boundaries_equal(neighbor_cluster, cluster))
          || (neighbor_cluster->boundary_vertexes[0] == contracted_v ||
          neighbor_cluster->boundary_vertexes[1] == contracted_v)){
          adjacencies[i] = cluster;
          vertex_added = true;
          break;
        }

      } 
    }
    if(!vertex_added && null_idx == -1){
      throw std::invalid_argument("Vertex could not be added to the tree");
    }
    if(!vertex_added){
      adjacencies[null_idx] = cluster;
    }
  }
}

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
  //tree and checking to see if atleast one neighbor of every vertex is in
  //the tree.
  for(auto vertex : *affected){
    if(contraction_tree[vertex][round] != nullptr){ 
      bool vertex_valid = false;
      for(int i = 0; i < degree_bound; i++){
        if(contraction_tree[vertex][round][i] != nullptr && 
          contraction_tree[vertex][round][i]->boundary_vertexes.size() == 2){
          auto dereferenced_cluster = *contraction_tree[vertex][round][i];
          int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
          if(maximal_set.count(vertex) == 0 && 
            ((contracts(neighbor, round, 0)  && affected->count(neighbor) == 0)|| 
            maximal_set.count(neighbor) == 1)){
            vertex_valid = true; // If either the vertex is next to an contracting unaffected neighbor or 
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
      if(maximal_set.count(vertex) == 0 && vertex_valid == false){
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
  bool spreads = true;
  for(int i = 0; i < degree_bound; i++){ //Neighbors of the neighbor of the affected vertex.
    auto neighbor_cluster = contraction_tree[vertex][round][i]; 
    if(neighbor_cluster != nullptr){
      RCCluster<aug_t> dereferenced_cluster = *neighbor_cluster;
      int neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster); // Neighbor of the affected vertex.
      // If my neighbor does not contract
      if(!(contracts(neighbor, round, 0) && current_affected->count(neighbor))) {
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
        if(cluster != nullptr){
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

  RCCluster<aug_t>* new_edge = (RCCluster<aug_t>*) new RCCluster<aug_t>(weight);
  new_edge->boundary_vertexes.push_back(u);
  new_edge->boundary_vertexes.push_back(v);

  add_neighbor(0, new_edge, -1); //v can be > t. Make it ternary.

  affected->insert(u);          // Insert initial affected vertices.
  affected->insert(v);


  representative_clusters[u] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[u]), PROCESSING};
  representative_clusters[v] = std::tuple<RCCluster<aug_t>*, int>{std::get<0>(representative_clusters[v]), PROCESSING};

  update(); 
}

template<typename aug_t>
void RCTree<aug_t>::cut(int u, int v) {
  //Delete an edge between 2 vertices of an RC Tree.
  //If no edge (u,v) exists, throw an exception.
}
/*
template<typename aug_t>
int RCTree<aug_t>::find_update_neighbor(vertex_t neighbor, vertex_t v, int round){
  Once the cluster of a contracting vertex determined, all neighbors to which it belongs need
 * to be updated with its new representative cluster. The method below finds the spot
 * in the neighbors adjacency list, where the new representative cluster belongs.

  auto adjacencies = contraction_tree[neighbor][round];
  int to_return = -1;
  for(int i = 0; i < degree_bound; i++){
    auto neighbor_ptr = contraction_tree[neighbor][round][i];

    // We check to see if the adjacency list has an empty spot
    // or if the rc cluster associated at that spot is from a
    // vertex that contracted earlier. If so, replace it
    // with the new RC Cluster being contracted.
    if(neighbor_ptr == nullptr ||
      representative_clusters[v] == neighbor_ptr || 
      (neighbor_ptr->round_contracted != PROCESSING 
      && neighbor_ptr->round_contracted < round)){
      to_return = i;
    }
  }
  if(to_return == -1){
    throw std::invalid_argument("Replacement index was not found");
  }
  return to_return;
}
*/

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
    new_cluster->aug_val += neighboring_cluster.aug_val; // Needs to be a associative binary op passed in by the user.
    if(neighboring_cluster.boundary_vertexes.size() == 2){
      neighbor = GET_NEIGHBOR(vertex, neighboring_cluster);
    }
  }

  new_cluster->boundary_vertexes.push_back(neighbor);
  new_cluster->parent = neighbor;
  //Add to the neighbor this new unary cluster.
  add_neighbor(round + 1, new_cluster, vertex);
  representative_clusters[vertex] = std::tuple{new_cluster, round}; //need to free the old tuple

  // Add other vertex I am raking onto into affected set if 
  // the cluster of the vertex being raked in the next round
  // is not a unary cluster. 
  // QUESTION : DO I NEED TO ADD IT IF VALUES NOT EQUIVALENT AS WELL - YES, to recontract with new values.

  affected->insert(neighbor);
  contraction_tree[vertex].resize(round + 1);
  representative_clusters[neighbor] = std::tuple{std::get<0>(representative_clusters[neighbor]), PROCESSING};
  affected->erase(vertex);
}

template<typename aug_t>
void RCTree<aug_t>::compress(vertex_t vertex, int round){

}

template<typename aug_t>
void RCTree<aug_t>::finalize(vertex_t v, int round){
  representative_clusters[v] = *(new std::tuple<RCCluster<aug_t>*, int>{(new RCCluster(0, -1)), round});
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
    //Insert new neighbor list for next round for each vertex. But this is inserting for those where
    //it is not even necessary and updating the list at the next round unnecessarily.
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


    for(auto vertex : *maximal_set){
      //Redo contractions for all vertices in maximal set.
      //Remove vertices from set of affected vertices.

      auto vertex_degree = get_degree(vertex, round);
      if(vertex_degree == 1) {
        rake(vertex, round);
      } else if(vertex_degree == 0) {
        // if the degree of a vertex is not 1 or 2 then it must be 0.
        // We have found the new root of the RC Tree and as such, we can
        // now finalize the vertex as the root of the new tree.
        finalize(vertex, round);
        break;
      }
    }
    
    std::unordered_set<int> new_affected;
    for(auto vertex : *affected){
      for(int i = 0; i < degree_bound; i++){
        auto cluster = contraction_tree[vertex][round][i];
        if(cluster != nullptr && cluster->boundary_vertexes.size() == 2) {
          auto dereferenced_cluster = *cluster;
          auto neighbor = GET_NEIGHBOR(vertex, dereferenced_cluster);
          /*  
          if(vertex == 2 && neighbor == 1 && round == 0){
            bool eq = std::get<0>(representative_clusters[vertex]) == contraction_tree[1][1][1];
            std::cout << std::get<0>(representative_clusters[vertex]) << std::endl;
            std::cout << contraction_tree[1][1][1] << std::endl;
            std::cout << eq << std::endl;
          }*/
          if(contracts(neighbor, round, 1)) continue;

          add_neighbor(round + 1, cluster, vertex);
          
          new_affected.insert(neighbor);
          std::get<1>(representative_clusters[neighbor]) = PROCESSING;
        }
      }
    }

    affected->merge(new_affected);
    round += 1;
    delete maximal_set;
  }
}
