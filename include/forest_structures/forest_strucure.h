#pragma once

#include "types.h"
#include "util.h"
#include <unordered_set>
#include <unordered_map>
#include <parlay/sequence.h>

using namespace std;
using namespace parlay;

class ForestStructure {
public:
    virtual void insert_vertices(sequence<vertex_t>& V) = 0;
    virtual void delete_vertices(sequence<vertex_t>& V) = 0;
    virtual void insert_edges(sequence<Edge>& E) = 0;
    virtual void delete_edges(sequence<Edge>& E) = 0;
    virtual sequence<bool> check_edges(sequence<Edge>& E) = 0;
    virtual vertex_t get_degree(vertex_t v) = 0;

    virtual sequence<vertex_t> get_parents(sequence<vertex_t>& V) = 0;
    virtual sequence<pair<vertex_t,vertex_t>> count_parents(sequence<vertex_t>& V) = 0;
    virtual void set_parents (sequence<vertex_t>& V, sequence<vertex_t>& P) = 0;
    virtual void unset_parents (sequence<vertex_t>& V) = 0;
    virtual void add_children(sequence<vertex_t>& V) = 0;
    virtual vertex_t get_child_count(vertex_t v) = 0;
};

class SequentialForest : public ForestStructure {
private:
    struct ForestNode {
        std::unordered_set<vertex_t> neighbors;
        vertex_t parent = -1;
        vertex_t child_count = 0;
        void insert_neighbor(vertex_t neighbor) { neighbors.insert(neighbor); }
        void remove_neighbor(vertex_t neighbor) { neighbors.erase(neighbor); }
        bool contains_neighbor(vertex_t neighbor) { return neighbors.find(neighbor) != neighbors.end(); }
    };

    std::unordered_map<vertex_t, ForestNode*> vertices;

public:
    void insert_vertices(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            auto node = new ForestNode();
            vertices[V[i]] = node;
        }
    }

    void delete_vertices(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            for (vertex_t neighbor : vertices[V[i]]->neighbors) {
                vertices[neighbor]->remove_neighbor(V[i]);
                vertices[V[i]]->remove_neighbor(neighbor);
            }
            delete vertices[V[i]];
            vertices.erase(V[i]);
        }
    }

    void insert_edges(sequence<Edge>& E) {
        for (int i = 0; i < E.size(); ++i) {
            vertices[E[i].src]->insert_neighbor(E[i].dst);
            vertices[E[i].dst]->insert_neighbor(E[i].src);
        }
    }

    void delete_edges(sequence<Edge>& E) {
        for (int i = 0; i < E.size(); ++i) {
            vertices[E[i].src]->remove_neighbor(E[i].dst);
            vertices[E[i].dst]->remove_neighbor(E[i].src);
        }
    }

    sequence<bool> check_edges(sequence<Edge>& E) {
        sequence<bool> output(E.size());
        for (int i = 0; i < E.size(); ++i) {
            if (vertices[E[i].src]->contains_neighbor(E[i].dst)) {
                assert(vertices[E[i].dst]->contains_neighbor(E[i].src));
                output[i] = true;
            } else {
                assert(!vertices[E[i].dst]->contains_neighbor(E[i].src));
                output[i] = false;
            }
        }
        return output;
    }
    vertex_t get_degree(vertex_t v) {
        return vertices[v]->neighbors.size();
    }

    sequence<vertex_t> get_parents(sequence<vertex_t>& V) {
        std::unordered_set<vertex_t> parents;
        for (int i = 0; i < V.size(); ++i) {
            parents.insert(vertices[V[i]]->parent);
        }
        auto output = sequence<vertex_t>(parents.size());
        int index = 0;
        for (auto parent : parents) {
            output[index++] = parent;
        }
        return output;
    }

    sequence<pair<vertex_t,vertex_t>> count_parents(sequence<vertex_t>& V) {
        std::unordered_map<vertex_t,vertex_t> parent_counts;
        for (int i = 0; i < V.size(); ++i) {
            if (parent_counts.find(vertices[V[i]]->parent) == parent_counts.end()) {
                parent_counts[vertices[V[i]]->parent] = 1;
            } else {
                parent_counts[vertices[V[i]]->parent] += 1;
            }
        }
        auto output = sequence<pair<vertex_t,vertex_t>>(parent_counts.size());
        int index = 0;
        for (auto entry : parent_counts) {
            output[index++] = entry;
        }
        return output;

    }

    void set_parents (sequence<vertex_t>& V, sequence<vertex_t>& P) {
        for (int i = 0; i < V.size(); ++i) {
            vertices[V[i]]->parent = P[i];
        }
    }

    void unset_parents (sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            vertices[V[i]]->parent = -1;
        }
    }

    void add_children(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            vertices[V[i]]->child_count += 1;
        }
    }

    vertex_t get_child_count(vertex_t v) {
        return vertices[v]->child_count;
    }
};
