#include <parlay>
#include<vector>

typedef vertex uint32_t;

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
    std::unordered_set<int> affected_vertices; 
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

int RCTree::get_degree(int u, int round) {
    
    int degree = 0;
    
    for(auto cluster : adj[round][u])
        if(cluster.boundaries.size() == 2) degree++;

    return degree;
}
RCTree::find_affected(int round) {
    
    std::queue<RCCluster> q(affected);
    std::unordered_set<int> visited(affected);

    while(!q.empty()) {
        
        RCCluster curr  = q.pop();
        
        if(visited.contains(curr)) {
            continue;
        }
        parlay::sequence<parlay::sequence<RCCluster>> induced_tree = adj[round];    
        
        for(int i = 0; i < 3; i++){
            
            // Found path cluster
            vector<int> boundaries = induced_tree[curr][i].boundary_vertices;
            if(boundaries.size() == 2) {
                // Find opposing vertex
                int neighbor = boundaries[0] != curr ? boundaries[0] : boundaries[1];
                
                if(induced_tree[neighbor].size() > 2) {
                    
                         
                }
            }
        }
    }
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
    
    affected_vertices.insert(u);          // Insert initial affected vertices.
    affected_vertices.insert(v);          // Better to use parlay::hashtable??
    
     
    update(); 
    
}

RCTree::cut(int u, int v) {
    //Delete an edge between 2 vertices of an RC Tree.
    //If no edge (u,v) exists, throw an exception.
    
    
}
RCTree::rake(int u) {

}

RCTree::compress(int u) {

}

RCTree::update() {
   // Based on set of affected vertices, find an MIS of the affected vertices, which induce
   // a subtree on the original tree and recontract them accordingly. Then determine new
   // set of affected vertices and recurse.
  

   // Sequential BFS running from each inital affected vertex at most distance 2.  
   BFS(affected_vertices, new_affected);
}




