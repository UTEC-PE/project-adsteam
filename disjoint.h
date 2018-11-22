//
// Created by Sebastian on 2/10/2018.
//

#ifndef GRAPH_DISJOINT_H
#define GRAPH_DISJOINT_H
#include <iostream>
#include <map>

using namespace std;
struct NodeDisjoint {
    char rank;
    char data;
    NodeDisjoint * parent;

    NodeDisjoint(char data) : data(data), rank(0), parent(this) {};
};

class DisjointSet {
private:
    map<char, NodeDisjoint *> nodes;

public:

    DisjointSet() {};

    void makeSet(char data) {
        NodeDisjoint *node = new NodeDisjoint(data);
        this->nodes[data] = node;
    }

    bool unionSet(char data1, char data2) {
        NodeDisjoint *parent1 = findSet(data1);
        NodeDisjoint *parent2 = findSet(data2);

        if (parent1 != parent2) {
            if (parent1->rank >= parent2->rank) {
                parent1->rank = (parent1->rank == parent2->rank) ? parent1->rank + 1 : parent1->rank;
                parent2->parent = parent1;
            } else {
                parent1->parent = parent2;
            }

            return true;
        }

        return false;
    }

    NodeDisjoint *findSet(char data) {
        return findSet(this->nodes[data]);
    }

    NodeDisjoint *findSet(NodeDisjoint *node) {
        NodeDisjoint *current = node;
        while (current != current->parent) {
            current = current->parent;
        }

        node->parent = current;
        return current;
    }

    char getParent(char data){
        return (nodes[data]->parent)->data;
    }

    ~DisjointSet() {}
};
#endif //GRAPH_DISJOINT_H