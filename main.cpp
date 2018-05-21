#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <memory>
#include "data.h"
#include "algo.h"
#include <unordered_map>
#include <omp.h>
//#include "point.h"

using namespace std;

//typedef std::unordered_map<int, double> Mymap;
//typedef std::unordered_map<pointD*, Mymap> Mymap1;
int main(int argc, char* args[]) {
    setlocale(LC_ALL, "Russian");
    try {
//        Mymap temp;
//        temp.insert(Mymap::value_type(1,5.6));
//        Mymap temp0;
//        temp0.insert(Mymap::value_type(5,7));
//        Mymap1 temp1;
//        pointD p(2.3,0.7);
//        pointD p2(2.3,0.7);
//        temp1.insert(Mymap1::value_type((&p),temp));
//        bool ex = false;
//        for(auto ent1 : temp1) {
//            if(*(ent1.first) == p2) {
//                ex = true;
//                temp1[ent1.first] = temp0;
//                break;
//            }
//        }
//        int a = 1;

//#pragma omp parallel for
//for(int i=0; i<10; i++) {
//    cout << i << endl;
//}
        algo::OpenFile(string(args[1]));
        try {
            algo::MainAlgorithm();
        }
        catch (string &exp) {
            cout << "Error: " << exp << endl;
        }
    }
    catch (string &exp) {
        cout << "An error occurred while attempting to open the file. The error is: " << exp << endl;
    }
    return 0;
}