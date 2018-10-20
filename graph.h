//
// Created by Sebastian on 28/09/2018.
//

#ifndef GRAPH_GRAPH_H
#define GRAPH_GRAPH_H
#include <iostream>
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include <string>
#include <set>
#include <stack>
#include <limits>
#include <algorithm>
#include "disjoint.h"

using namespace std;

struct Edge
{
    char start;
    char final;
    int weight;
    void printEdge()
    {
        cout << "{" << start << "," << final << "}" << endl;
    }

};


class Graph
{
private:
    typedef map<char, vector<Edge>> adjList;
    adjList graphmap;
    int vertices;
    int edges;
    bool dir;


public:

    Graph(bool direccionado=true)
    {
        vertices= 0;
        edges=0;
        dir = direccionado;

    };

    Graph(fstream& Document)
    {
        edges = 0;
        vertices = 0;

        vector<string> datum;
        string word;
        word.clear();

        while (Document >> word)
        {
            datum.push_back(word);
        }

        int v = stoi(datum[0]);
        dir = stoi(datum[1]);

        int n = stoi(datum[2]);

        for(int k = 0; k < v; k++)
        {
            insertNode(datum[k + 3].c_str()[0]);
        }

        int j = v + 3;

        for(int i = 0; i < n; i++)
        {
            addEdge(datum[j].c_str()[0], datum[j+1].c_str()[0], stoi(datum[j+2]));
            j += 3;
        }
    }

    void printGraph(){
        for (auto pair: graphmap){
            cout << pair.first << ": ";
            for (auto edge: pair.second){
                cout << " ";
                edge.printEdge();
                cout << " ";
            }
            cout << endl;
        }
    }

    void insertNode(char node){
        vector<Edge> edges;
        auto pair= make_pair(node, edges);
        graphmap.insert(pair);
        vertices++;
    };

    void addEdge(char start, char final, int weight) {
        Edge edge{start, final, weight};

        if (graphmap.count( start ) && graphmap.count( final )) {
            graphmap[start].push_back(edge);
            if(!dir) {
                Edge edge2{final, start, weight};
                graphmap[final].push_back(edge2);
                edges++;
            }
            edges++;
        }
        else{
            cout << "Primero debe insertar los nodos";
        }
    };

    void removeEdge(char start, char final){

        int counter=0;
        for (auto edge: graphmap[start]){
            if(edge.final==final) graphmap[start].erase(graphmap[start].begin()+counter);
            counter++;
        }
        edges--;
        if(!dir){
            for (auto edge: graphmap[final]){
                if(edge.final==start) graphmap[final].erase(graphmap[final].begin()+counter);
                counter++;
            }
            edges--;
        }

    };

    void removeNode(char start){
        if (dir){
            graphmap.erase(start);
            vertices--;
            for (auto pair: graphmap){
                auto node= pair.first;
                removeEdge(node, start);
            }
        }
        else{
            for (auto edge: graphmap[start]){
                auto connectedNode= edge.final;
                removeEdge(connectedNode, start);
            }
            graphmap.erase(start);
            vertices--;
        }

    };

    bool findVertex(char vertex){
        for (auto pair : graphmap){
            if (pair.first== vertex){
                return true;
            }
        }
        return false;

    };
    bool findEdge(char start, char final){
        if (findVertex(start) && findVertex(final)){
            for (auto edge: graphmap[start]){
                if (edge.final== final) return true;
            }
        }
        return false;

    };

    Edge findMinimumEdge(vector<char> &VisitedVertex)
    {
        int min = numeric_limits<int>::max();
        Edge MinEdge;
        for(int i = 0; i < VisitedVertex.size(); i++)
        {
            for (int j = 0; j < graphmap[VisitedVertex[i]].size(); j++) {
                if((graphmap[VisitedVertex[i]][j].weight < min) && (count(VisitedVertex.begin(), VisitedVertex.end(), graphmap[VisitedVertex[i]][j].final ) == 0))
                {
                    MinEdge = graphmap[VisitedVertex[i]][j];
                    min = MinEdge.weight;
                }
            }
        }
        VisitedVertex.push_back(MinEdge.final);
        return MinEdge;
    }


    void primAlgorithm(char StartingPoint)
    {
        vector<char> VisitedVertex;
        vector<Edge> MST;
        VisitedVertex.push_back(StartingPoint);
        while(VisitedVertex.size() != vertices)
        {
            MST.push_back(findMinimumEdge(VisitedVertex));
        }

        cout << "{";
        for (auto &i : MST) {
            i.printEdge();
        }
        cout << "}";
    }

    Edge findSmallestEdge(vector<char> &VisitedVertex)
    {
        int min = numeric_limits<int>::max(); //Infinity
        Edge MinEdge;

        for (auto f: graphmap)
        {
            for(auto j: f . second)
            {
                if(j.weight < min && count(VisitedVertex.begin(), VisitedVertex.end(), j.final) == 0)
                {
                    min = j.weight;
                    MinEdge = j;
                }
            }
        }
        if (VisitedVertex.size() == 0) VisitedVertex.push_back(MinEdge.start);
        VisitedVertex.push_back(MinEdge.final);
        return MinEdge;
    }


    void kruskalAlgorithm()
    {
        vector<char> VisitedVertex;
        vector<Edge> MST;
        while(VisitedVertex.size() != vertices)
            MST.push_back(findSmallestEdge(VisitedVertex));
        cout << "{" << endl;
        for (auto &i : MST) {
            i.printEdge();
        }
        cout << "}";
    }

    bool isConnected(){
        DisjointSet disjointSet;

        for (auto pair: graphmap){
            auto vertex= pair.first;
            disjointSet.makeSet(vertex);
        }

        for (auto pair: graphmap){
            for (auto edge: pair.second){
                disjointSet.unionSet(edge.start, edge.final);
            }
        }
        set<char> parents;
        for (auto pair: graphmap){
            auto vertex= pair.first;
            parents.insert(disjointSet.getParent(vertex));
        }
        if(parents.size()==1){
            cout << "Todos estan enlazados";
            return true;
        }else{
            cout << "No están conectados";
            return false;
        }

    }


    void depthSons(char start, set<char> &visitedNodes, stack<char> &stackClosedNodes){
        visitedNodes.insert(start);
        for (auto edge: graphmap[start]){
            if (visitedNodes.find(edge.final)==visitedNodes.end()){
                depthSons(edge.final, visitedNodes, stackClosedNodes);
            }
        }
        stackClosedNodes.push(start);
    };


    bool isStronglyConnected(){
        if (!dir){
            cout << "Esta operacion solo es para grafos dirigidos" << endl;
            return false;
        }

        if (isConnected()){
            char startNode= (*(graphmap.begin())).first;
            set<char> visitedNodes;
            stack<char> stackClosedNodes;
            for (auto pair: graphmap) {
                auto currentNode= pair.first;
                if(visitedNodes.find(currentNode) ==visitedNodes.end()) {
                    depthSons(startNode, visitedNodes, stackClosedNodes);
                }
            }

            vector<char> closedNodes;
            while(!stackClosedNodes.empty()){
                closedNodes.push_back(stackClosedNodes.top());
                stackClosedNodes.pop();
            }
            DisjointSet disjointSet;

            for (auto node: closedNodes){
                disjointSet.makeSet(node);
            }

            for (auto vertex: closedNodes){
                for (auto edge: graphmap[vertex]){
                    disjointSet.unionSet(edge.final, edge.start);
                }
            }

            set<char> parents;
            for (auto node: closedNodes){
                parents.insert(disjointSet.getParent(node));
            }

            if(parents.size()==1){
                cout << "fuertemente conexo";
                return true;
            }else{
                cout << "no es fuertemente conexo";
                return false;
            }

        }
        return false;
    };




    bool isBipartite(){
        //grafos no conexos?
        map<char, bool> nodeAndColor;
        for (auto pair: graphmap){
            auto newPair= make_pair(pair.first, NULL);
            nodeAndColor.insert(newPair);
        }
        nodeAndColor.begin()->second=true;
        auto start=  nodeAndColor.begin()->first;
        vector<char> bfsNodes;
        set<char> visitedNodes;
        bfsNodes.push_back(start);

        while(!bfsNodes.empty()){
            auto current=*bfsNodes.begin();
            for (auto edge: graphmap[current] ){
                auto nextNode= edge.final;
                if(nodeAndColor[nextNode]==NULL){
                    nodeAndColor[nextNode]= !nodeAndColor[current];
                }
                else{
                    if(nodeAndColor[nextNode]==nodeAndColor[current]) return false;
                }
                if (visitedNodes.find(nextNode) == visitedNodes.end())
                {
                    bfsNodes.push_back(nextNode);
                    visitedNodes.insert(nextNode);
                }
            }
            bfsNodes.erase(bfsNodes.begin());
        }
        return true;
    }

    bool depth(char currentNode, char final, set<char>& visitedNodes){

        bool flag= false;
        visitedNodes.insert(currentNode);
        if (currentNode==final) return true;
        for (auto edge: graphmap[currentNode]){
            auto nextNode= edge.final;
            if (visitedNodes.find(nextNode) == visitedNodes.end()) {
                flag= depth(nextNode, final, visitedNodes);
            }
            if (flag) break;
        }
        return flag;
    };

    bool DFS(char start, char final){
        set<char> visitedNodes;
        return depth(start, final, visitedNodes);
    };


    bool BFS(char start, char final){
        vector<char> bfsNodes;
        set<char> visitedNodes;
        bfsNodes.push_back(start);
        visitedNodes.insert(start);

        while(!bfsNodes.empty()){
            auto current= *bfsNodes.begin();
            for (auto edge: graphmap[current]){
                auto nextNode= edge.final;

                if(nextNode==final) return true;

                if (visitedNodes.find(nextNode) == visitedNodes.end())
                {
                    bfsNodes.push_back(nextNode);
                    visitedNodes.insert(nextNode);
                }
            }
            bfsNodes.erase(bfsNodes.begin());
        }
        return false;
    };


    float density(){
        return (float(edges)/float((vertices*(vertices-1))));
    }
    bool isDense(float x= 0.7){ //un valor por defecto
        return density()>x;
    }
    bool isSparse(float x= 0.7){ //un valor por defecto
        return density()<x;
    }

    int vertexGrade(char start){
        return graphmap[start].size();
    }

    bool isSource(char start){
        if(dir) {
            for (auto f: graphmap) {
                for (auto i: f.second) {
                    if (i.final == start) return false;
                }
            }
            return true;
        }
        else{
            cout << "Esta operacion es para grafos direccionados";
            return false;
        }
    }

    bool isLeaf(char start){
        if(dir) {
            return graphmap[start].empty();
        }
        else{
            cout << "Esta operacion es para grafos direccionados";
            return false;
        }
    }


    ~Graph(){
        vertices=0;
        edges=0;
        graphmap.clear();
    }

};

#endif //GRAPH_GRAPH_H
