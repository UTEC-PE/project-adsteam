//
// Created by Sebastian on 2/10/2018.
//

#ifndef GRAPH_DISJOINT_H
#define GRAPH_DISJOINT_H
#include <iostream>
#include <map>

using namespace std;
struct Node {
    char rank;
    char data;
    Node* parent;

    Node(char data) : data(data), rank(0), parent(this) {};
};

class DisjointSet {
private:
    map<char, Node *> nodes;

public:

    DisjointSet() {};

    void makeSet(char data) {
        Node *node = new Node(data);
        this->nodes[data] = node;
    }

    bool unionSet(char data1, char data2) {
        Node *parent1 = findSet(data1);
        Node *parent2 = findSet(data2);

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

    Node *findSet(char data) {
        return findSet(this->nodes[data]);
    }

    Node *findSet(Node *node) {
        Node *current = node;
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
