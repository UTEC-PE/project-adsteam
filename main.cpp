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

    cout << "Test de floydWarshall" << endl;
    cout << "---------------------" << endl;

    TestGraph.floydWarshall();

    cout << "Test de bellmanFord" << endl;
    cout << "-------------------" << endl;

    TestGraph.bellmanFord(&Zero);  

    cout << endl;

    cout << "Test de Djikstra" << endl;
    cout << "----------------" << endl;

    TestGraph.djikstraAlgorithm(&Zero, &One);

    cout << endl;

    cout << "Test de A*" << endl;
    cout << "----------" << endl;

    TestGraph.aStarAlgorithm(&Zero,&One);

    cout << endl;

    cout << "Test de GreedyBFS" << endl;
    cout << "-----------------" << endl;

    TestGraph.greedyBFS(&Zero, &Two);

}