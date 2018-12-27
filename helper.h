//
// Created by ANDGRA on 31.05.2018.
//

#ifndef HELPER_H
#define HELPER_H

#include <iostream>
#include <math.h>
#include <algorithm>
#include <string>
#include <sstream>  //for std::istringstream
#include <vector>  //for std::istringstream
#include <typeinfo>
#include <unordered_map>
//#include "comp.h"
#include <iomanip>
#include "simple_fft/fft_settings.h"
#include "simple_fft/fft.h"

using std::count;
using std::string;
using std::getline;
using std::vector;
using std::stringstream;
using std::to_string;
using std::unordered_map;
using std::locale;
using std::setprecision;
using std::fixed;

typedef complex_type comp;

string left_trim(const string& str, char symbol = ' ')
{
    size_t first = str.find_first_not_of(symbol);
    if (string::npos == first)
    {
        return str;
    }
    return str.substr(first);
}

string right_trim(const string& str, char symbol = ' ')
{
    size_t first = str.find_first_not_of(symbol);
    if (string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(symbol);
    if (last != first) {
        last = (last - first + 1);
    }
    return str.substr(0, last + 1);
}

string trim(const string& str, char symbol = ' ')
{
    return left_trim(right_trim(str, symbol), symbol);
}

vector<string> split(const string s, char delim)
{
    vector<string> elems;
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim))
    {
        elems.push_back(item);
    }
    return elems;
}

constexpr unsigned int str2int(const char* str, int h = 0)
{
    return !str[h] ? 5381 : (str2int(str, h+1) * 33) ^ str[h];
}
template<typename K, typename V>

bool containsKey(unordered_map<K, V> m, K k, V& v) {
    auto it = m.find(k);
    if(it != m.end()) {
        v = it->second;
        return true;
    }
    return false;
}

template<typename K, typename V>
bool containsKey(unordered_map<K, V> m, K k) {
    auto it = m.find(k);
    return it != m.end();
}

class My_punct : public std::numpunct<char> {
protected:
    char do_decimal_point() const {return ',';}//comma
};

double stod_c(string s) {
    locale loc(locale(), new My_punct);
    stringstream ss(s);
    ss.imbue(loc);
    double d;
    ss >> d;
    return d;
}

template<typename T>
string serialize(T* arr, int size) {
    string res = "[";
    for (int i = 0; i < size; i++) {
        res += to_string(arr[i]);
        if (i < size - 1) {
            res += ",";
        }
    }
    res += "]";
    return res;
}


string to_string(comp numb) {
    stringstream ss1, ss2;
    string re, im;
    ss1 << fixed << setprecision( 17 ) << numb.real();
    ss2 << fixed << setprecision( 17 ) << numb.imag();
    ss1 >> re;
    ss2 >> im;
    return "{" + re + "," + im + "}";
}


void conjugate(comp *array, unsigned int size)
{
    for(int i = 0; i < size; i++)
        array[i] = conj(array[i]);
}

void FFT(comp* data, unsigned int size) {
    const char * error = NULL; // error description
    simple_fft::FFT(data, size, error);
//    conjugate(data, size);
}

void IFFT(comp* data, unsigned int size) {
    const char * error = NULL; // error description
    simple_fft::IFFT(data, size, error, false);
}

#if defined(_WIN32)
#include <iostream>
#include <chrono>

class Timer
{
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return std::chrono::duration_cast<second_>
                (clock_::now() - beg_).count(); }

private:
    typedef std::chrono::high_resolution_clock clock_;
    typedef std::chrono::duration<double, std::ratio<1> > second_;
    std::chrono::time_point<clock_> beg_;
};

#else
#include <iostream>
#include <ctime>

class Timer
{
public:
    Timer() { clock_gettime(CLOCK_REALTIME, &beg_); }

    double elapsed() {
        clock_gettime(CLOCK_REALTIME, &end_);
        return end_.tv_sec - beg_.tv_sec +
               (end_.tv_nsec - beg_.tv_nsec) / 1000000000.;
    }

    void reset() { clock_gettime(CLOCK_REALTIME, &beg_); }

private:
    timespec beg_, end_;
};
#endif

#endif //HELPER_H
