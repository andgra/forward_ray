#include <iostream>
#include <string>
#include "algo.h"

using namespace std;

int main(int argc, char* args[]) {
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
    return 0;
}