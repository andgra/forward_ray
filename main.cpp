#include <iostream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <memory>
//#include "data.h"
//#include "algo.h"
#include <unordered_map>
#include "point.h"

using namespace std;

typedef std::unordered_map<int, double> Mymap;
typedef std::unordered_map<pointD*, Mymap> Mymap1;
int main(int argc, char* args[]) {
    try {
        Mymap temp;
        temp.insert(Mymap::value_type(1,5.6));
        Mymap1 temp1;
        pointD p(2.3,0.7);
        pointD p2(2.3,0.7);
        temp1.insert(Mymap1::value_type((&p),temp));
        bool ex = false;
        for(auto ent1 : temp1) {
            if(*(ent1.first) == p2) {
                ex = true;
                break;
            }
        }
        int a = 1;
//        algo::OpenFile(string(args[1]));
        try {
//            data.MainAlgorithm();
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