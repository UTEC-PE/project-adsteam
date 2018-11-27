#include <iostream>
#include <fstream>
#include <cstdlib>
#include "graphProject.h"
#include <vector>
#include <map>

int main() 
{

    cout << "Prueba de Algoritmos" << endl;
    cout << "--------------------" << endl;

    cout << endl;

    cout << "Creacion de Nodos" << endl;
    cout << "------------------" << endl;

    Node<int> Zero(0,0,2);
    Node<int> One(1,1,2);
    Node<int> Two(2,2,2);
    Node<int> Three(3,3,2);
    Node<int> Four(4,4,1);
    Node<int> Five(5,3,0);
    Node<int> Six(6,2,0);
    Node<int> Seven(7,1,0);
    Node<int> Eight(8,2,1);

    cout << endl;

    cout << "Inserción de Nodos" << endl;
    cout << "------------------" << endl;

    Graph<int> TestGraph(false);

    TestGraph.insertNode(&Zero);
    TestGraph.insertNode(&One);
    TestGraph.insertNode(&Two);
    TestGraph.insertNode(&Three);
    TestGraph.insertNode(&Four);
    TestGraph.insertNode(&Five);
    TestGraph.insertNode(&Six);
    TestGraph.insertNode(&Seven);
    TestGraph.insertNode(&Eight);

    TestGraph.printGraph();

    cout << endl;

    cout << "Insercion de Edges" << endl;
    cout << "------------------" << endl;

    TestGraph.addEdge(&Zero,&One,4);
    TestGraph.addEdge(&Zero,&Seven,8);
    TestGraph.addEdge(&One,&Two,8);
    TestGraph.addEdge(&One,&Seven,11);
    TestGraph.addEdge(&Two,&Three,7);
    TestGraph.addEdge(&Two,&Eight,2);
    TestGraph.addEdge(&Two,&Five,4);
    TestGraph.addEdge(&Three,&Five,14);
    TestGraph.addEdge(&Three,&Four,9);
    TestGraph.addEdge(&Four,&Five,10);
    TestGraph.addEdge(&Five,&Six,2);
    TestGraph.addEdge(&Six,&Eight,6);
    TestGraph.addEdge(&Six,&Seven,1);
    TestGraph.addEdge(&Seven,&Eight,7);

    cout << endl;

    TestGraph.printGraph();

    TestGraph.bellmanFord(&Zero);


    /*    Node<char> A('A',1,3);
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

        cout << "Test de Heuristica" << endl;
        cout << "------------------" << endl;
        cout << "La Distancia entre los nodos A y B es: " << TestGraph.aStarHeuristic(&A, &D);

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