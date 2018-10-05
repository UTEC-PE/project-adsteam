#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Graph.h"
#include <vector>
#include <map>
using namespace std;


int main() {


    vector<int> v;
    v.push_back(10);
    v.push_back(10);
    v.push_back(10);
    v.push_back(10);
    v.push_back(30);
    int i=0;
    for (auto f: v){
        cout << f;
    }
    cout << endl;
    for (auto f: v){
        if (f==10) {v.erase(v.begin()+i);} i++;}
        i=0;
    for (auto f: v){
        if (f==10) {v.erase(v.begin()+i);} i++;}
        i=0;
    for (auto f: v){
        if (f==10) {v.erase(v.begin()+i);} i++;}
    for (auto f: v){
        cout << f;
    }


        return 0;
}