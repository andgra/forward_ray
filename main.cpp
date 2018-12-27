#include <iostream>
#include <string>
#include "algo.h"
#include "point.h"
#include "helper.h"
#include <unordered_map>
//#include <Exocor>
#include "mpi.h"

using std::cout;
using std::string;
using std::endl;

typedef unordered_map<int, string> mm;

int main(int argc, char* args[]) {
//    Exocortex.DSP.Complex a(123);
//    pointI p(100, 30);
//
//    auto p1 = pointF::conv(p);
//    int a = 1;
//    auto s = pointD::serialize(p);
//    auto a = pointD::deserialize(s);
//    auto s1 = pointD::serialize(a);
//    auto b = s == s1;
//    auto b1 = p == a;
//    mm m;
//    m.insert(mm::value_type(1,"asd"));
////    auto it = m.find(2);
//    string* s;
//    for(auto it: m) {
//        s = &(m[it.first]);
//    }
//    *s = "fsd";
//    delete s;
//    auto c = containsKey(m, 1, s);

    MPI_Init(&argc, &args);
    try {
        if(args[1] == nullptr) {
            throw string("parameter was not sent");
        }
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

    MPI_Finalize();
    return 0;
}