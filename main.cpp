#include <iostream>
#include <string>
//#include "algo.h"

using namespace std;

int main(int argc, char* args[]) {
    auto i = (int)1;
    cout << to_string(i) << endl;
    try {
//        algo::OpenFile(string(args[1]));
        try {
//            algo::MainAlgorithm();
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