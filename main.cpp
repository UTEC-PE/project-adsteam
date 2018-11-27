#include <iostream>
#include <fstream>
#include <cstdlib>
#include "graphProject.h"
#include <vector>
#include <map>

int main() 
{

        Node<char> A('A',1,2);
        Node<char> B('B',2,2);
        Node<char> D('D',3,3);

        Edge<char> C(&A,&B,3);

        Graph<char> TestGraph;
        TestGraph.insertNode(&A);
        TestGraph.insertNode(&B);
        TestGraph.insertNode(&D);

        TestGraph.addEdge(&A,&B,4);
        TestGraph.removeEdge(&B,&A);
        TestGraph.addEdge(&A,&D,7);

        TestGraph.removeNode(&B);

        TestGraph.printGraph();

        cout << "Test de Find:" << endl;
        if(TestGraph.findVertex(&A))
            cout << "Se encontro A." << endl;
        if(!(TestGraph.findVertex(&B)))
            cout << "No se encontro B." << endl;

        cout << "Test de Find Edge:" << endl;
        if(TestGraph.findEdge(&A,&D))
            cout << "Se encontro el Edge AD" << endl;
        if(!TestGraph.findEdge(&D,&B))
            cout << "No se encontro el Edge DB" << endl;

        TestGraph.primAlgorithm(&A);
        TestGraph.kruskalAlgorithm();

        if(TestGraph.isBipartite())
            cout << "Es Bipartito." << endl;
        else
            cout << "No es Bipartito" << endl;

        /*graph.insertNode('A',1,1);
        graph.insertNode('B',2,2);
        graph.insertNode('C',3,3);
        graph.insertNode('D',4,4);
        graph.addEdge('A', 'B', 6);
        graph.addEdge('C', 'D', 1);
        graph.addEdge('C', 'B', 1);
        graph.addEdge('D', 'B', 5);*/

        //fstream Document;
        //Document.open("graphStart.txt");
        //Graph Test(Document);


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
}