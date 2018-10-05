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
#include "disjoint.h"

using namespace std;

/*El archivo que lee este grafo sigue el formato:
 * 11 1
 * 1 2 3
 * 1 3 5
 * O sea, hay 11 vertices, es True (1) que sea direccionado, el 1 conecta con 2 con un peso 3, con 3 con un peso 5.
*/

struct Edge
{
    char start;
    char final;
    int weight;

};


//asume que el grafo es direccionado

class Graph
{
private:
    typedef map<char, vector<Edge>> adjList;
    adjList graphmap;
    int vertices;
    int edges;
    bool dir;

    bool depth(char start, char final, set<char>& my_set){

        bool i= false; //si no encuentra nada, deberá devolver falso
        my_set.insert(start); //inserta el nodo start para no volver a visitarlo
        if (start==final) return true; //si estás en el que buscas, retornas true
        for (auto f: graphmap[start]){ //f es cada nodo conectado a tu start
            auto it = my_set.find(f.final); //busca si el nodo conectado ya ha sido visitado
            if (it == my_set.end()) i= depth(f.final, final, my_set); //si no ha sido visitado, visitar
            if (i) break; //si el nodo visitado devolvio true, deja de buscar y retorna true.
        } //si se acaban los hijos o estas en un ciclo, deja de moverte y retorna el valor i por defecto (false).
        return i;
    };

public:

    Graph()
    {
        vertices= 0;
        edges=0;
        dir= true;

    };

    Graph(fstream Document)
    {
        edges = 0;

        Document.open("graphStart.txt"); //Abrimos el archivo que va a tener este nombre
        vector<string> datum; //Aqui van a ir todas los numeros que diga el archivo.
        string word; //Cada uno de los datos
        word.clear();

        while (Document >> word)
        {
            datum.push_back(word); //Transforma los string en numero y los pone en el datum
        } //Mientras hayan palabras en el documento. Leera cada una y la pushback en el vector de datos.

        dir = stoi(datum[1]);
        vertices = stoi(datum[0]);
        int j = 3;

        for(int i = 0; i < vertices; i++)
        {
            addEdge(datum[j].c_str()[0], datum[j+1].c_str()[0], stoi(datum[j+2]));
            //There might be a more efficient way of doing this.
            j += 3;
        }
    }

    //Primero insertas nodos, luego los conectas. No es simultáneo.

    void insertNode(char start){

        vector<Edge> edges;
        auto pair= make_pair(start, edges);
        graphmap.insert(pair);
        vertices++;

    };

    void addEdge(char start, char final, int weight) {
        Edge edge{start, final, weight};

        //asume que es direccionado
        if (graphmap.count( start ) && graphmap.count( final )) {

            if (graphmap[start].empty()) {
                graphmap[start].push_back(edge);
            } else {
                graphmap[start].push_back(edge);
            }

            edges++;
        }
        else{
            cout << "Primero debe insertar los nodos";
        }
    };

    void removeEdge(char start, char final){
        int i=0;
        for (auto f: graphmap[start]){
            if(f.final==final) graphmap[start].erase(graphmap[start].begin()+i);
            i++;
        }
        edges--;
    };


    void removeNode(char start){

        graphmap.erase(start);
        vertices--;

        for (auto f: graphmap){
        removeEdge(f.first, start);
        }



    };
    bool is_connected(){
        //usar disjoint
        DisjointSet my_set;

        for (auto f: graphmap){
            my_set.makeSet(f.first); //iniciar cada vertice como un nodo
        }

        for (auto f: graphmap){
            for (auto i: f.second){ //unir cada edge, los que esten conectados apuntarán al mismo padre.
                my_set.unionSet(i.start, i.final);
            }
        }

        set check;
        for (auto f: graphmap){ //acumula los padres de cada nodo
            check.insert(my_set.get_Parent(f.first));
        }
        return check.size()==1; //si solo hay un padre, es porque todos están conectados.

    }

    bool not_connected(){
        return !is_connected();
    }

    //FALTA ESTO: ???
    bool strong_connected(){
        //primer enfoque: si hay conexión entre todos y A está conectado a todos, el camino es ida y vuelta.
        if(is_connected()){
            return (graphmap.begin()->second == vertices-1);
        }
        return false;
        //segundo enfoque: Usar BFS o DFS para saber si están conectados en ambos sentidos.
    }


    bool findVertex(char start){
        for (auto f : graphmap){ //si el nodo ha sido inicializado...
            if (f.first== start){
                return true;
            }
        }
        return false;

    };
    bool findEdge(char start, char final){
        if (findVertex(start)){ //si existe el nodo, busca
            for (auto f: graphmap[start]){
                if (f.final== final) return true;
            }
        }
        return false;

    };

    bool is_bipartite(){ //se puede dividir en dos conjuntos?
        map<char, bool> my_map; //vamos a pintar cada nodo con false o true
        for (auto f: graphmap){
            auto pair= make_pair(f.first, NULL); //cada pareja empieza con NULL
            my_map.insert(pair);
        }
        my_map.begin()->second=true; //empezamos con true
        for (auto f: graphmap){
            for (auto i:f.second){
                if (my_map[i.final]==NULL){ //si no ha sido visitado...
                    my_map[i.final]= !my_map[i.first]; //ponerle el contrario (distinto color)
                }
                else{ //si ha sido visitado y es igual al color del vecino, devuelve false
                    if(my_map[i.final]== my_map[i.first]){return false; }
                }
            }
        }
        //si cada uno se puede pintar de color distinto al vecino, devuelve true.
        return true;
    }



    bool DFS(char start, char final){
        set<char> my_set; //crear un set donde vas a poner los nodos que vayas visitando.
      return depth(start, final, my_set); //va a aplicar recursividad
    };


    bool BFS(char start, char final){
        set<char> my_set; //inicializa un set de nodos visitados.
        my_set.insert(start); //inserta el primer nodo
        while(!my_set.empty()){ //mientras queden nodos por visitar, avanza
            auto current= graphmap[*my_set.begin()]; //el primer nodo
            for (auto f: current ){ //para cada hijo del nodo
                if(f==final) return true; //si algún hijo es el final, retorna true
                my_set.insert(f); //los hijos van a la cola
            }
            my_set.erase(my_set.begin());  //borramos el nodo que ha sido visitado, seguimos.
        }
        return false; //si has revisado a todos los descendientes y no hay nada, retorna false.
    };


    float density(){
        //asumiendo dirigido
        return (edges/(vertices*(vertices-1)));
    }
    bool isDense(float x= 0.7){ //un valor por defecto
        return density()>x;
    }
    bool isDispersed(float x= 0.7){ //un valor por defecto
        return density()<x;
    }

    int vertex_grade(char start){
        return graphmap[start].size();
    }

    bool is_root(char start){  //solo nacen nodos de él?

        //busca por todos los vectores. Si start no aparece como final ni una vez, retorna true
        for (auto f: graphmap){
            for (auto i: f.second){
                if(i.final == start) return false;
            }
        }
        return true;
    }

    bool is_leaf(char start){ //solo hay nodos hacia él?
        return graphmap[start].empty(); //si no tiene nigun saliente, asumimos que sí.
    }


    ~Graph(){
        vertices=0;
        edges=0;
        graphmap.clear();
    }

};

#endif //GRAPH_GRAPH_H
