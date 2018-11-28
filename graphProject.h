#ifndef GRAPHPROJECT_H
#define GRAPHPROJECT_H
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
#include <math.h>
#include "disjoint.h"

using namespace std;

template <typename T>

struct Node
{
    T name;
    int xCoordinate;
    int yCoordinate;

    Node(T nameSetup, int xSetup, int ySetup)
    {
        name = nameSetup;
        xCoordinate = xSetup;
        yCoordinate = ySetup;
    }

    void relocate(int newXCoordinate, int newYCoordinate)
    {
        xCoordinate = newXCoordinate;
        yCoordinate = newYCoordinate;
    }
};

template <typename T>

struct Edge
{
    Node<T> * start;
    Node<T> * final;
    int weight;

    Edge()
    {
        start = nullptr;
        final = nullptr;
        weight = 0;
    }

    Edge(Node<T> * startingNode, Node<T> * endingNode, int peso)
    {
        start = startingNode;
        final = endingNode;
        weight = peso;
    }

    void printEdge()
    {
        cout << "{" << start -> name << "," << final -> name << "}" << endl;
    }
};

template <typename T>

struct Path
{
    vector<Edge<T>> edgePath;
    int weight = 0;

    Path()
    {
        weight = numeric_limits<int>::max();
    }

    Path(int weightSetup)
    {
        weight = weightSetup;
    }

    void setVector(vector<Edge<T>> vectorSetup)
    {
        edgePath = vectorSetup;
    }

    void setWeight(int value)
    {
        weight = value;
    }

    void calculateWeight()
    {
        weight = 0;

        for(int i = 0; i < edgePath.size(); i++)
            weight += edgePath[i].weight;
    };

    void addEdge(Edge<T> edgeToPush)
    {
        edgePath.push_back(edgeToPush);
        calculateWeight();
    }

    void reversal()
    {
        reverse(edgePath.begin(),edgePath.end());
    }

    void printPath()
    {
        for (auto f: edgePath)
            f.printEdge();
    }

};

template <typename T>

class Graph
{
private:
    typedef map<Node<T> *, vector<Edge<T>>> adjList;
    adjList graphmap;
    int vertices;
    int edges;
    bool dir;

public:

    Graph(bool direccionado = false)
    {
        vertices = 0;
        edges = 0;
        dir = direccionado;
    };

    /*Graph(fstream& Document)
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
    }*/

    void printGraph()
    {
        if (!vertices) 
            cout << "El grafo no tiene vertices T_T.";
        for (auto pair : graphmap)
        {
            cout << (pair.first) -> name << ": ";
            for (auto edge : pair.second)
            {
                cout << " ";
                edge.printEdge();
                cout << " ";
            }
            cout << endl;
        }
    };

  void insertNode(Node<T> * nodeToInsert)
    {
        vector<Edge<T>> edges;
        auto pair = make_pair(nodeToInsert, edges);
        if(graphmap.find(nodeToInsert) == graphmap.end())
            vertices++;
        graphmap.insert(pair);
    };

    void addEdge(Node<T> * start, Node<T> * final, int weight)
     {
        Edge<T> edgeToAdd{start, final, weight};

        if (graphmap.count(start) && graphmap.count(final)) 
        {
            graphmap[start].push_back(edgeToAdd);
            if(!dir) 
            {
                Edge<T> edgeToAdd2{final, start, weight};
                graphmap[final].push_back(edgeToAdd2);
                edges++;
            }
            edges++;
        }
        else
        {
            cout << "Los nodos puestos como parametros no estan incluidos en la lista de adyacencia.";
        }
    };

    void removeEdge(Node<T> * start, Node<T> * final)
    {
        int counter=0;
        for (auto edge : graphmap[start]){
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

    void removeNode(Node<T> * start)
    {
        if (dir)
        {
            graphmap.erase(start);
            vertices--;
            for (auto pair: graphmap){
                auto node = pair.first;
                removeEdge(node, start);
            }
        }

        else
        {
            for (auto edge: graphmap[start]){
                auto connectedNode = edge.final;
                removeEdge(connectedNode, start);
            }
            graphmap.erase(start);
            vertices--;
        }
    };

    bool findVertex(Node<T> * vertex)
    {
        for (auto pair : graphmap){
            if (pair.first == vertex){
                return true;
            }
        }
        return false;
    };

    bool findEdge(Node<T> * start, Node<T> * final)
    {
        if (findVertex(start) && findVertex(final))
        {
            for (auto edge: graphmap[start])
            {
                if (edge.final == final) return true;
            }
        }
        return false;
    };

    Edge<T> findMinimumEdge(vector<Node<T> *> &VisitedVertex)
    {
        int min = numeric_limits<int>::max();
        Edge<T> MinEdge;
        for(int i = 0; i < VisitedVertex.size(); i++)
        {
            for (int j = 0; j < graphmap[VisitedVertex[i]].size(); j++) 
            {
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


    void primAlgorithm(Node<T> * startingPoint)
    {
        vector<Node<T> *> VisitedVertex;
        vector<Edge<T>> MST;
        VisitedVertex.push_back(startingPoint);
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

    Edge<T> findSmallestEdge(vector<Node<T> *> &VisitedVertex)
    {
        int min = numeric_limits<int>::max();
        Edge<T> MinEdge;

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
        vector<Node<T> *> VisitedVertex;
        vector<Edge<T>> MST;
        while(VisitedVertex.size() != vertices)
            MST.push_back(findSmallestEdge(VisitedVertex));
        cout << "{" << endl;
        for (auto &i : MST) {
            i.printEdge();
        }
        cout << "}";
    }

    bool isConnected()
    {
        DisjointSet disjointSet;

        for (auto pair: graphmap){
            auto vertex = (pair.first)->name;
            disjointSet.makeSet(vertex);
        }

        for (auto pair: graphmap){
            for (auto edge: pair.second){
                disjointSet.unionSet((edge.start)->name, (edge.final)->name);
            }
        }
        set<T> parents;
        for (auto pair: graphmap){
            auto vertex = (pair.first)->name;
            parents.insert(disjointSet.getParent(vertex));
        }
        if(parents.size()==1){
            cout << "Es un grafo conexo.";
            return true;
        }else{
            cout << "No es un grafo conexo.";
            return false;
        }
}


    void depthSons(T start, set<T> &visitedNodes, stack<T> &stackClosedNodes)
    {
        visitedNodes.insert(start);
        for (auto edge: graphmap[start]){
            if (visitedNodes.find((edge.final)->name)==visitedNodes.end()){
                depthSons((edge.final)->name, visitedNodes, stackClosedNodes);
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
            T startNode= ((*(graphmap.begin())).first)->name;
            set<T> visitedNodes;
            stack<T> stackClosedNodes;
            for (auto pair: graphmap) {
                auto currentNode= (pair.first)->name;
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

    bool isBipartite()
    {
        map<Node<T> *, bool> nodeAndColor;

        for (auto pair: graphmap)
        {
            auto newPair= make_pair(pair.first, NULL);
            nodeAndColor.insert(newPair);
        }

        nodeAndColor.begin() -> second = true;
        auto start =  nodeAndColor.begin()->first;
        vector<Node<T> *> bfsNodes;
        set<Node<T> *> visitedNodes;
        bfsNodes.push_back(start);

        while(!bfsNodes.empty())
        {
            auto current = *bfsNodes.begin();
            for (auto edge: graphmap[current] ){
                auto nextNode= edge.final;
                if(!(nodeAndColor[nextNode])){
                    nodeAndColor[nextNode] = !nodeAndColor[current];
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
    };

    bool depth(Node<T> * currentNode, Node<T> * final, set<Node<T> *>& visitedNodes, vector<Edge<T>>& tree)
    {

        bool flag = false;
        visitedNodes.insert(currentNode);
        if (currentNode==final) 
            return true;
        for (auto edge: graphmap[currentNode])
        {
            auto nextNode= edge.final;
            if (visitedNodes.find(nextNode) == visitedNodes.end()) 
            {
                Edge<T> edge{currentNode, nextNode, edge.weight};
                tree.push_back(edge);
                flag= depth(nextNode, final, visitedNodes, tree);
            }
            if (flag) break;
        }
        return flag;
    };

    bool DFS(Node<T> * start, Node<T> * final, vector<Edge<T>>& tree)
    {
        set<Node<T> *> visitedNodes;
        return depth(start, final, visitedNodes, tree);
    };

    void print_DFS(Node<T> * start, Node<T> * final)
    {
        vector<Edge<T>> tree;
        DFS(start, final, tree);
        for (auto f: tree)
        {
            f.printEdge();
            cout << endl;
        }
    };


    bool BFS(Node<T> * start, Node<T> * final, vector<Edge<T>>& tree){
        vector<Node<T> *> bfsNodes;
        set<Node<T> *> visitedNodes;
        bfsNodes.push_back(start);
        visitedNodes.insert(start);
        while(!bfsNodes.empty()){
            auto current= *bfsNodes.begin();
            for (auto edge: graphmap[current]){
                auto nextNode= edge.final;
                Edge<T> *tedge = new Edge<T>{current, nextNode, edge.weight};
                tree.push_back(*tedge);
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

    void print_BFS(Node<T> * start, Node<T> * final){
        vector<Edge<T>> tree;
        BFS(start, final, tree);
        for (auto f: tree){
            f.printEdge();
            cout << endl;
        }
    }


    float density(){
        return (float(edges)/float((vertices*(vertices-1))));
    }
    bool isDense(float x= 0.7){ //un valor por defecto
        return density()>x;
    }
    bool isSparse(float x= 0.7){ //un valor por defecto
        return density()<x;
    }

    int vertexGrade(Node<T> * start){
        return graphmap[start].size();
    }

    bool isSource(Node<T> * start){
        if(dir) 
        {
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

    bool isLeaf(Node<T> * start){
        if(dir) 
        {
            return graphmap[start].empty();
        }
        else{
            cout << "Esta operacion es para grafos direccionados";
            return false;
        }
    }

    //Segunda Entrega de Projecto.
    Node<T>* nameToNode(T name){
        for(auto f: graphmap){
            if (f.first->name == name) {
                auto node= f.first;
                return node;
            }
        }
    }

    Edge<T> findEdgeWith(Node<T> * nodePrev, vector<Edge<T>> MST)
    {
        for(int i = 0; i < MST.size(); i++)
        {
            if(MST[i].final == nodePrev)
            {
                return MST[i];
            }
        }
    }

    Path<T> * greedyBFS(Node<T> * sourceNode, Node<T>* finalNode)
    {
        vector<Node<T> *> VisitedVertex;
        vector<Edge<T>> MST;
        VisitedVertex.push_back(sourceNode);

        do
        {
            MST.push_back(findMinimumEdge(VisitedVertex));
        }while(MST[MST.size() - 1].final != finalNode);

        Path<T> * greedyPath = new Path<T>();
        Edge<T> edgeToAdd;

        edgeToAdd = MST[MST.size()-1];
        greedyPath -> addEdge(edgeToAdd);

        while(edgeToAdd.start != sourceNode)
        {
            greedyPath -> addEdge(edgeToAdd);
            edgeToAdd = findEdgeWith(edgeToAdd.start, MST);
            cout << edgeToAdd.start -> name << " " << edgeToAdd.final -> name << endl;
        }

        greedyPath -> reversal();

        return greedyPath;
    }


    map<T, int> bellmanFord(Node<T>* Node){

        map<T, int> table;

        set<T> reach;

        reach.insert(Node->name);

        bool flag=true;
        int n;
        while(flag){
            n=reach.size();
            for (auto i: reach){

                for(auto j: graphmap[nameToNode(i)]){
                    reach.insert((j.final)->name);
                }
            }
            if(n==reach.size()){
                flag=false;
            }
        }

        auto firstPair= make_pair(Node->name, 0);
        table.insert(firstPair);
        int inf= numeric_limits<int>::max();

        for (auto pair: graphmap){
            if(pair.first->name!= Node->name && reach.count((pair.first)->name)){
                auto newPair= make_pair((pair.first)->name, inf);
                table.insert(newPair);
            }
        }

        for (int i=0; i<vertices-1; i++){
            for (auto vertex: table){
                for(auto neighbor: graphmap[nameToNode(vertex.first)]){
                    if((vertex.second)!=inf){
                        if (vertex.second + neighbor.weight< table[(neighbor.final)->name]){
                            table[(neighbor.final)->name]= vertex.second + neighbor.weight;
                        }
                    }
                }
            }
        }
        cout << endl << "Matriz: " << endl;
        for (auto f: table){
            cout << f.first << ": " << f.second << endl;
        }
        return table;
    }


    pair<map<Node<T>*, map<Node<T>*, int>>, map<Node<T>*, map<Node<T>*, Node<T>*>>>  floydWarshall() {
        map<Node<T>*, map<Node<T>*, int>> distanceMatrix;
        map<Node<T>*, map<Node<T>*, Node<T>*>> pathMatrix;
        int inf;
        for (auto vertice:  graphmap){
            map<Node<T>*, int> vertice2;
            map<Node<T>*, Node<T>*> vertice2b;

            auto pair= make_pair(vertice.first, vertice2);
            auto pairb= make_pair(vertice.first, vertice2b);
            distanceMatrix.insert(pair);
            pathMatrix.insert(pairb);


            for (auto verticeB: graphmap){
                map<Node<T>*, int> temp;
                inf = numeric_limits<int>::max();
                auto pair2 = make_pair(verticeB.first, inf);
                auto pair2b= make_pair(verticeB.first, vertice.first);
                distanceMatrix[vertice.first].insert(pair2);
                pathMatrix[vertice.first].insert(pair2b);
            }

            distanceMatrix[vertice.first][vertice.first]= 0;
            pathMatrix[vertice.first][vertice.first]= vertice.first;

        }

        for (auto vertice: graphmap){
            for (auto neighbor: vertice.second){
                distanceMatrix[vertice.first][neighbor.final]= neighbor.weight;
            }
        }


        for (auto verticeI: graphmap){
            for(auto verticeJ: graphmap){
                for (auto verticeK: graphmap) {
                    if (inf != distanceMatrix[verticeJ.first][verticeI.first] &&
                        inf != distanceMatrix[verticeI.first][verticeK.first]) {
                        auto sum = distanceMatrix[verticeJ.first][verticeI.first] +
                                   distanceMatrix[verticeI.first][verticeK.first];
                        if (sum < distanceMatrix[verticeJ.first][verticeK.first]) {
                            distanceMatrix[verticeJ.first][verticeK.first] = sum;
                            pathMatrix[verticeJ.first][verticeK.first] = verticeI.first;
                        }

                    }
                }
            }
        }

        cout << "Distance: " << endl;
        for (auto i: distanceMatrix){
            cout << (i.first)->name << "| ";
            for (auto j: i.second){
                cout << (j.first)->name << ": "<< j.second << " ";
            }
            cout << endl;
        }
        cout << endl;
        cout <<"Path: " << endl;
        for (auto i: pathMatrix){
            cout << (i.first)->name << "| ";
            for (auto j: i.second){
                cout << (j.first)->name << ":"<< (j.second)->name<<" ";
            }
            cout << endl;
        }
        auto my2Tables= make_pair(distanceMatrix, pathMatrix);
        return my2Tables;
    }

    void includeNode(vector<Node<T>*> & visitedNodes,  vector<Node<T>*> & nonVisitedNodes, Node<T> * nodeToInclude)
    {
        nonVisitedNodes.erase(remove(nonVisitedNodes.begin(), nonVisitedNodes.end(), nodeToInclude), nonVisitedNodes.end());
        visitedNodes.push_back(nodeToInclude);
    }

    bool checkIfLess(vector<Node<T>*> & visitedNodes, Node<T> * nodeToCheck, map<Node<T> *, Path<T> *> shortestValueTables, int counter)
    {
        return (!(count(visitedNodes.begin(), visitedNodes.end(), graphmap[nodeToCheck][counter].final)) //If not in the visited nodes.
        && (graphmap[nodeToCheck][counter].weight + (shortestValueTables[nodeToCheck] -> weight) < shortestValueTables[graphmap[nodeToCheck][counter].final] -> weight));
        }

    void updateTable(map<Node<T> *, Path<T> *> &shortestValueTables, Edge<T> edgeToAdd)
    {   
        shortestValueTables[edgeToAdd.final] -> setWeight(shortestValueTables[edgeToAdd.start] -> weight);
        shortestValueTables[edgeToAdd.final] -> setVector(shortestValueTables[edgeToAdd.start] -> edgePath);
        shortestValueTables[edgeToAdd.final] -> addEdge(edgeToAdd);
    }

    void djikstraSmallestPath(vector<Node<T> *> visitedNodes, Node<T> * smallestNodeFound, map<Node<T> *, Path<T>*> &shortestValueTables)
    {
        for(int i = 0; i < graphmap[smallestNodeFound].size(); i++)
        {
            if(checkIfLess(visitedNodes, smallestNodeFound, shortestValueTables, i))
            {
                updateTable(shortestValueTables, graphmap[smallestNodeFound][i]);
            }
        }

    }

    Node<T> * calculateMinimum(map<Node<T> *, Path<T> *> shortestValueTables, vector<Node<T> *> nonVisitedNodes)
    {
        int min = numeric_limits<int>::max();
        Node<T> * minNode = nullptr;
        for(int i = 0; i < nonVisitedNodes.size(); i++)
        {
            //cout << shortestValueTables[nonVisitedNodes[i]] -> weight << endl;
            if(shortestValueTables[nonVisitedNodes[i]] -> weight < min)
            {
                min = shortestValueTables[nonVisitedNodes[i]] -> weight;
                minNode = nonVisitedNodes[i];
            }
        }
        return minNode;
    }

    Path<T> * djikstraAlgorithm(Node<T> * sourceNode, Node<T> * finalNode)
    {
        map<Node<T> *, Path<T> *> shortestValueTables;
        vector<Node<T> *> visitedNodes;
        vector<Node<T> *> nonVisitedNodes;

        Node<T> * smallestNodeFound = sourceNode;

        for(auto f: graphmap)
        {
            shortestValueTables[f.first] = new Path<T>();
            nonVisitedNodes.push_back(f.first);
        }

        shortestValueTables[sourceNode] -> setWeight(0);
        includeNode(visitedNodes, nonVisitedNodes, sourceNode);
        
        do
        {
            djikstraSmallestPath(visitedNodes, smallestNodeFound, shortestValueTables);
            smallestNodeFound = calculateMinimum(shortestValueTables, nonVisitedNodes);
            includeNode(visitedNodes, nonVisitedNodes, smallestNodeFound);
        }
        while(smallestNodeFound);

        for(auto node: shortestValueTables)
        {
            cout << node.first -> name << ": ";
            cout << node.second -> weight << endl;
        }

        return shortestValueTables[finalNode];
    }

    double aStarHeuristic(Node<T> * sourceNode, Node<T> * finalNode)
    {
        return sqrt(pow((finalNode -> yCoordinate) - (sourceNode -> yCoordinate), 2) + pow((finalNode -> xCoordinate) - (sourceNode -> xCoordinate), 2));

    }

    bool aStarcheckIfLess(vector<Node<T>*> & visitedNodes, Node<T> * nodeToCheck, map<Node<T> *, Path<T> *> shortestValueTables, int counter)
    {
        return (!(count(visitedNodes.begin(), visitedNodes.end(), graphmap[nodeToCheck][counter].final))
         && (graphmap[nodeToCheck][counter].weight + (shortestValueTables[nodeToCheck]) -> weight) < shortestValueTables[graphmap[nodeToCheck][counter].final]);
    }

    void aStarSmallestPath(vector<Node<T> *> visitedNodes, Node<T> * smallestNodeFound, map<Node<T> *, Path<T> *> & shortestValueTables)
    {
        for(int i = 0; i < graphmap[smallestNodeFound].size(); i++)
        {
            if(checkIfLess(visitedNodes, smallestNodeFound, shortestValueTables, i))
                updateTable(shortestValueTables, graphmap[smallestNodeFound][i]);
        }
    }

    Node<T> * aStarCalculateMinimum(map<Node<T> *, Path<T> *> shortestValueTables, vector<Node<T> *> nonVisitedNodes, Node<T> * finalNode)
    {
        int min = numeric_limits<int>::max();
        Node<T> * minNode = nullptr;
        for(int i = 0; i < nonVisitedNodes.size(); ++i)
        {
            if(shortestValueTables[nonVisitedNodes[i]] -> weight + aStarHeuristic(nonVisitedNodes[i], finalNode) < min) //Vas por buen camino continua con esto.
            {
                min = shortestValueTables[nonVisitedNodes[i]] -> weight + aStarHeuristic(nonVisitedNodes[i], finalNode);
                minNode = nonVisitedNodes[i];
            }
        }

        return minNode;
    }

    Path<T> * aStarAlgorithm(Node<T> * sourceNode, Node<T> * finalNode)
    {
        map<Node<T> *, Path<T> *> shortestValueTables;
        vector<Node<T> *> visitedNodes;
        vector<Node<T> *> nonVisitedNodes;

        Node<T> * smallestNodeFound = sourceNode;

        for(auto f: graphmap)
        {
            shortestValueTables[f.first] = new Path<T>();
            nonVisitedNodes.push_back(f.first);
        }

        shortestValueTables[sourceNode] -> setWeight(0);
        includeNode(visitedNodes, nonVisitedNodes, sourceNode);
        
        do
        {
            aStarSmallestPath(visitedNodes, smallestNodeFound, shortestValueTables);
            smallestNodeFound = aStarCalculateMinimum(shortestValueTables, nonVisitedNodes, finalNode);
            includeNode(visitedNodes, nonVisitedNodes, smallestNodeFound);
        }
        while(smallestNodeFound != finalNode || !(smallestNodeFound));

        for(auto node: shortestValueTables)
        {
            cout << node.first -> name << ": ";
            cout << node.second -> weight << endl;
        }

        return shortestValueTables[finalNode];
    }

    ~Graph()
    {
        vertices = 0;
        edges = 0;
        graphmap.clear();
    }

};

#endif //GRAPHPROJECT_H