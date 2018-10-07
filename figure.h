//
// Created by ANDGRA on 20.05.2018.
//

#ifndef FIGURE_H
#define FIGURE_H

#include <vector>
#include "edge.h"
#include <limits>
#include "point.h"

using std::vector;
using std::numeric_limits;
using std::min;
using std::max;

class figure//фигура задаётся набором отрезков
{
public:
    int index;
    vector<edge> edges;//координаты фигур целочисленны
    int speed;
    double density;
    double absorption_c;
    double gamma;
    bool isItDifObject(){ return _isItDifObjest; }
    int left;
    int right;
    int bottom;
    int up;
    double area;

    vector<pointI> points()
    {
        auto res = vector<pointI>(edges.size());
        for(edge e: edges) {
            res.push_back(e.p1);
        }
        return /*Array.ConvertAll<Point, PointF>(*/res; //, (p=>(PointF)p));
    }

    figure() {}

    figure(vector<edge> edges_i, int speed_i, double density_i, double absorption_i, int index_i) {
        edges = edges_i;
        speed = speed_i;
        density = density_i;
        absorption_c = absorption_i;
        index = index_i;
        gamma = density * speed;
        _isItDifObjest = (GetArea() < 36) && (index != 0);
        left = numeric_limits<int>::max();
        bottom = numeric_limits<int>::max();
        right = numeric_limits<int>::min();
        up = numeric_limits<int>::min();
        for(edge e: edges)
        {
            int curLeft = min(e.p1.X, e.p2.X);
            int curBottom = min(e.p1.Y, e.p2.Y);
            int curRight = max(e.p1.X, e.p2.X);
            int curTop = max(e.p1.Y, e.p2.Y);
            if(curLeft < left) left = curLeft;
            if(curBottom < bottom) bottom = curBottom;
            if(curRight > right) right = curRight;
            if(curTop > up) up = curTop;
        }
        area = GetArea();
    }

    bool ContainsPoint(pointD point) { //bottleneck
        bool result = false;
        if (index == 0)
            return true;
        if (point.X > right || point.X < left || point.Y > up || point.Y < bottom)
            return false;

        int cnt = edges.size();
        int j = cnt - 1;
        for (int i = 0; i < cnt; i++) {
            if (edges[i].p2.Y < point.Y && edges[j].p2.Y >= point.Y ||
                edges[j].p2.Y < point.Y && edges[i].p2.Y >= point.Y) {
                if (edges[i].p2.X +
                    (point.Y - edges[i].p2.Y) / (edges[j].p2.Y - edges[i].p2.Y) * (edges[j].p2.X - edges[i].p2.X) <=
                    point.X) {
                    result = !result;
                }
            }
            j = i;
            if (edges[i].ContainsPoint(point))//тот алгоритм не обрабатывает часть граничных случаев ( пара 50 - 75)
            {
                //succeeded++;
                return true;
            }
        }
        //if (result == true)
        //    succeeded++;
        //else
        //    failed++;
        return result;
    }

    double radius = 0.5;

    bool IsPointDif(pointD point) {
        for(edge e: edges)
        {
            if (GetDistance(point, pointD::conv(e.p1)) < radius)
                return true;
        }
        return false;
    }

    pointD GetNearestDifPoint(pointD point) {
        double minDist = numeric_limits<double>::max();
        double curDist;
        int index = -1;
        for (int i = 0; i < edges.size(); i++) {
            curDist = GetDistance(point, pointD::conv(edges[i].p1));
            if (minDist > curDist) {
                minDist = curDist;
                index = i;
            }
        }
        return pointD::conv(edges[index].p1);
    }

private:
    bool _isItDifObjest;

    double GetDistance(pointD p1, pointD p2) {
        double difX = p2.X - p1.X;
        double difY = p2.Y - p1.Y;
        return sqrt(difX * difX + difY * difY);
    }

    double GetArea()//площадь полигона
    {
        double sum = 0;
        for(edge e: edges)
        {
            sum += double((e.p2.X - e.p1.X) * (e.p2.Y + e.p1.Y)) / 2;
        }
        return fabs(sum);
    }
};

#endif //FIGURE_H
