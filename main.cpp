#include <stdio.h>
#include <conio.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <memory>

using namespace std;

int main() {try
    {
        //MyData::OpenFile(args[1]);
        try
        {
            //MyData::MainAlgorithm();
            throw string("534");
        }
        catch (string& exp)
        {
            cout << "Error: " << exp << endl;
        }
    }
    catch (string& exp)
    {
        cout << "An error occurred while attempting to open the file. The error is: " << exp << endl;
    }
    return 0;
}