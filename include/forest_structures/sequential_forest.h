#pragma once

#include "bridge.h"
#include <forest_strucure.h>
#include <unordered_set>
#include <unordered_map>


using namespace parlay;

class SequentialNeighborIterator : public NeighborIterator {
private:
    std::unordered_set<vertex_t>::iterator iter;
    std::unordered_set<vertex_t>::iterator end;
public:
    SequentialNeighborIterator(std::unordered_set<vertex_t>::iterator _iter, std::unordered_set<vertex_t>::iterator _end) : iter(_iter), end(_end) {}
    vertex_t next() {
        if (iter == end) {
            return NONE;
        } else {
            return *(iter++);
        }
    }
};

class SequentialForest : public ForestStructure {
private:
    struct ForestNode {
        std::unordered_set<vertex_t> neighbors;
        vertex_t child_count = 0;
        vertex_t parent = NONE;
        vertex_t partner = NONE;
        uint32_t priority;
        ClusterStatus status = NORMAL;
        void insert_neighbor(vertex_t neighbor) { neighbors.insert(neighbor); }
        void remove_neighbor(vertex_t neighbor) { neighbors.erase(neighbor); }
        bool contains_neighbor(vertex_t neighbor) { return neighbors.find(neighbor) != neighbors.end(); }
    };

    std::unordered_map<vertex_t, ForestNode*> vertices;

public:
    void insert_vertices(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            if (V[i] == NONE) continue;
            if (vertices.find(V[i]) != vertices.end()) continue;
            auto node = new ForestNode();
            node->priority = hash32(V[i]);
            vertices[V[i]] = node;
        }
    }

    void delete_vertices(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            for (vertex_t neighbor : vertices[V[i]]->neighbors) {
                vertices[neighbor]->remove_neighbor(V[i]);
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

    sequence<Edge> map_edges_to_parents(sequence<Edge>& E) {
        sequence<Edge> output;
        for (int i = 0; i < E.size(); ++i) {
            if (vertices[E[i].src]->parent != NONE && vertices[E[i].dst]->parent != NONE) {
                if (vertices[E[i].src]->parent != vertices[E[i].dst]->parent) {
                    output.push_back({vertices[E[i].src]->parent,vertices[E[i].dst]->parent});
                }
            }
        }
        return output;
    }

    sequence<vertex_t> get_parents(sequence<vertex_t>& V) {
        std::unordered_set<vertex_t> parents;
        for (int i = 0; i < V.size(); ++i) {
            if (vertices[V[i]]->parent != NONE) {
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
            vertices[V[i]]->parent = NONE;
        }
    }

    void add_children(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            if (V[i] == NONE) continue;
            vertices[V[i]]->child_count += 1;
        }
    }

    void subtract_children(sequence<vertex_t>& V) {
        for (int i = 0; i < V.size(); ++i) {
            if (V[i] == NONE) continue;
            vertices[V[i]]->child_count -= 1;
        }
    }

    vertex_t get_degree(vertex_t v) {
        return vertices[v]->neighbors.size();
    }

    std::unique_ptr<NeighborIterator> get_neighbor_iterator(vertex_t v) {
        return std::make_unique<SequentialNeighborIterator>(vertices[v]->neighbors.begin(), vertices[v]->neighbors.end());
    }

    vertex_t get_first_neighbor(vertex_t v) {
        return *vertices[v]->neighbors.begin();
    }

    vertex_t get_other_neighbor(vertex_t v, vertex_t x) {
        auto iter = vertices[v]->neighbors.begin();
        vertex_t first = *iter;
        iter++;
        vertex_t second = (iter == vertices[v]->neighbors.end()) ? NONE : *iter;
        if (first != x) return first;
        return second;
    }

    vertex_t get_parent(vertex_t v) {
        return vertices[v]->parent;
    }

    void set_parent (vertex_t v, vertex_t p) {
        vertices[v]->parent = p;
    }

    void unset_parent(vertex_t v) {
        vertices[v]->parent = NONE;
    }

    vertex_t get_child_count(vertex_t v) {
        return vertices[v]->child_count;
    }

    bool contracts(vertex_t v) {
        if (vertices[v]->parent != NONE) {
            for (vertex_t neighbor : vertices[v]->neighbors) {
                if (vertices[neighbor]->parent == vertices[v]->parent) {
                    return true;
                }
            }
        }
        return false;
    }

    vertex_t get_partner(vertex_t v) {
        return vertices[v]->partner;
    }

    void set_partner(vertex_t v, vertex_t p) {
        vertices[v]->partner = p;
    }

    bool try_set_partner_atomic(vertex_t v, vertex_t p) {
        if (vertices[v]->partner == NONE) {
            return gbbs::CAS(&vertices[v]->partner, NONE, p);
        }
        return false;
    }

    void unset_partner(vertex_t v) {
        vertices[v]->partner = NONE;
    }

    bool is_local_max_priority (vertex_t v) {
        for (vertex_t neighbor : vertices[v]->neighbors) {
            if (vertices[neighbor]->parent != NONE || get_degree(neighbor) != 2) continue;
            if (vertices[neighbor]->priority >= vertices[v]->priority) {
                return false;
            }
        }
        return true;
    }

    ClusterStatus get_status(vertex_t v) {
        return vertices[v]->status;
    }

    void set_status(vertex_t v, ClusterStatus s) {
        vertices[v]->status = s;
    }

    bool try_set_status_atomic(vertex_t v, ClusterStatus s) {
        if (vertices[v]->status == NORMAL) {
            return gbbs::CAS(&vertices[v]->status, NORMAL, s);
        }
        return false;
    }

    void unset_status(vertex_t v) {
        vertices[v]->status = NORMAL;
    }

    void mark(vertex_t v) {
        set_partner(v, MARK);
    }

    void unmark(vertex_t v) {
        unset_partner(v);
    }

    bool is_marked(vertex_t v) {
        return vertices[v]->partner == MARK;
    }
};
