#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Graph.h"
#include <vector>
#include <map>
using namespace std;


int main() {

        Graph graph;
        graph.insertNode('A');
        graph.insertNode('B');
        graph.insertNode('C');
        graph.insertNode('D');
        //graph.insertNode('E');
        graph.addEdge('A', 'B', 6);
        graph.addEdge('B', 'A', 1);
        graph.addEdge('C', 'D', 1);
        graph.addEdge('C', 'B', 1);
        graph.addEdge('D', 'B', 5);
        Graph ciclo;
        ciclo.insertNode('A');
        ciclo.insertNode('B');
        ciclo.insertNode('C');
        ciclo.addEdge('A', 'B', 5);
        ciclo.addEdge('B', 'C', 5);
        ciclo.addEdge('C', 'A', 5);

        /*probando delete
        graph.removeEdge('B', 'A');
        graph.printGraph();
        graph.removeNode('B');
        graph.printGraph();
         */
        /*
        cout << graph.findVertex('B') << graph.findEdge('B', 'A');
        cout << graph.isConnected();
         cout << ciclo.isConnected()<< ciclo.isStronglyConnected() << ciclo.isBipartite();
         cout << ciclo.DFS('A', 'C') << ciclo.BFS('B', 'A') << graph.DFS('D', 'C');
         cout << ciclo.density() << graph.density();
         cout << graph.vertexGrade('C') << graph.isLeaf('C') << graph.isSource('C');
        */


        ciclo.primAlgorithm('A');
        cout << endl << endl;
        ciclo.kruskalAlgorithm();
        cout << endl << endl;
        graph.primAlgorithm('A');
        cout << endl << endl;
        graph.kruskalAlgorithm();





        return 0;
}