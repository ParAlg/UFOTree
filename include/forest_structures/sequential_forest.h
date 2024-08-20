#pragma once

#include <forest_strucure.h>
#include <unordered_set>
#include <unordered_map>


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

    sequence<vertex_t> get_endpoints(sequence<Edge>& E) {
        std::unordered_set<vertex_t> endpoints;
        for (int i = 0; i < E.size(); ++i) {
            endpoints.insert(E[i].src);
            endpoints.insert(E[i].dst);
        }
        auto output = sequence<vertex_t>(endpoints.size());
        int index = 0;
        for (auto endpoint : endpoints) {
            output[index++] = endpoint;
        }
        return output;
    }

    sequence<Edge> filter_edges(sequence<Edge>& E) {
        sequence<Edge> output;
        for (int i = 0; i < E.size(); ++i) {
            if (vertices[E[i].src]->parent != vertices[E[i].dst]->parent) {
                output.push_back(E[i]);
            }
        }
        return output;
    }

    sequence<vertex_t> get_parents(sequence<vertex_t>& V) {
        std::unordered_set<vertex_t> parents;
        for (int i = 0; i < V.size(); ++i) {
            if (vertices[V[i]]->parent != -1) {
                parents.insert(vertices[V[i]]->parent);
            }
        }
        auto output = sequence<vertex_t>(parents.size());
        int index = 0;
        for (auto parent : parents) {
            output[index++] = parent;
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
            if (V[i] == -1) continue;
            vertices[V[i]]->child_count += 1;
        }
    }

    void subtract_children(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            if (V[i] == -1) continue;
            vertices[V[i]]->child_count -= 1;
        }
    }

    vertex_t get_degree(vertex_t v) {
        return vertices[v]->neighbors.size();
    }

    sequence<vertex_t> get_neighbors(vertex_t v) {
        sequence<vertex_t> output;
        for (auto neighbor : vertices[v]->neighbors)
            output.push_back(neighbor);
        return output;
    }

    vertex_t get_parent(vertex_t v) {
        return vertices[v]->parent;
    }

    vertex_t get_child_count(vertex_t v) {
        return vertices[v]->child_count;
    }

    bool contracts(vertex_t v) {
        for (vertex_t neighbor : vertices[v]->neighbors) {
            if (vertices[neighbor]->parent == vertices[v]->parent) {
                return true;
            }
        }
        return false;
    }
};