//
// Created by ANDGRA on 20.05.2018.
//

#ifndef POINT_H
#define POINT_H
#include "helper.h"

template<class T>
class point
{
public:
    T X, Y;
    point() : X(0), Y(0) {}
    point(T x, T y) : X(x), Y(y) {}
    static point conv(point<int> p) {
        return point(p.X,p.Y);
    }
    static point conv(point<double> p) {
        return point(p.X,p.Y);
    }
    static point conv(point<float> p) {
        return point(p.X,p.Y);
    }
    bool operator==(point p) {
        return p.X == X && p.Y == Y;
    }
    bool operator < (const point & p_lhs) const
    {
        if(p_lhs.X < X) { return true ; }
        if(p_lhs.X > X) { return false ; }
        return (p_lhs.Y < Y) ;
    }
    bool operator > (const point & p_lhs) const
    {
        return !(this == p_lhs) && !(this < p_lhs);
    }

    static string serialize(point p) {
        return "" + to_string(p.X) + ";" + to_string(p.Y);
    }

    static point deserialize(string s) {
        vector<string> coords = split(s, ';');
        switch(str2int(typeid(T).name())) {
            case str2int("d"): return point(stod(coords[0]),stod(coords[1]));
            case str2int("f"): return point(stof(coords[0]),stof(coords[1]));
            case str2int("i"): return point(stoi(coords[0]),stoi(coords[1]));
            default: return point(0,0);
        }
    }
};

typedef point<float> pointF;
typedef point<double> pointD;
typedef point<int> pointI;

#endif //POINT_H
