//
// Created by ANDGRA on 20.05.2018.
//

#ifndef STRAIGHT_RAY_EDGE_H
#define STRAIGHT_RAY_EDGE_H

#include "vec2.h"
#include "point.h"

using std::min;
using std::max;

class edge {
public:
    pointI pStart;//начало отрезка
    pointI pEnd;//конец отрезка
    int edgeIndex;//номер грани
    int figureIndex;//номер фигуры
    int xStart, xEnd, yStart, yEnd;

    pointI p1;//на всякий случай, чтобы ничего не сломать. Разберусь - исправлю
    pointI p2;

    int A, B, C;//коэффициенты уравнения прямой

    vec2f vector() const{return vec2f::normalize(vec2f(p2.X - p1.X, p2.Y - p1.Y));}
    vec2f vectorOrdered() const{return vec2f::normalize(vec2f(xEnd - xStart, yEnd - yStart));}
    vec2f normalVector() const{return vec2f::normalize(vec2f(A, B));}

    edge(const pointI& _p1, const pointI& _p2, int eindex, int findex) {
        this->p1 = _p1;
        this->p2 = _p2;
        A = p1.Y - p2.Y;
        B = p2.X - p1.X;
        C = p1.X * p2.Y - p2.X * p1.Y;

        //направление слева направо(при перпендикуляре снизу вверх)
        if (p1.X < p2.X) {
            pStart = p1;
            pEnd = p2;
        } else if (p1.X > p2.X) {
            pStart = p2;
            pEnd = p1;
        } else {
            if (p1.Y < p2.Y) {
                pStart = p1;
                pEnd = p2;
            } else {
                pStart = p2;
                pEnd = p1;
            }
        }
        xStart = pStart.X;
        xEnd = pEnd.X;
        yStart = min(p1.Y, p2.Y);
        yEnd = max(p1.Y, p2.Y);
        this->edgeIndex = eindex;
        this->figureIndex = findex;
    }

    bool ContainsPoint(pointF p) {
        //if (p.X > xEnd || p.X < xStart || p.Y > yEnd || p.Y < yStart)
        //    return false;
        //float t = 1E-03f;
        //// ensure points are collinear
        //var zero = (p2.X - p1.X) * (p.Y - p1.Y) - (p.X - p1.X) * (p2.Y - p1.Y);
        //if (zero > t || zero < -t) return false;

        //// check if x-coordinates are not equal
        //if (p1.X - p2.X > t || p2.X - p1.X > t)
        //    // ensure x is between a.x & b.x (use tolerance)
        //    return p1.X > p2.X
        //        ? p.X + t > p2.X && p.X - t < p1.X
        //        : p.X + t > p1.X && p.X - t < p2.X;

        //// ensure y is between a.y & b.y (use tolerance)
        //return p1.Y > p2.Y
        //    ? p.Y + t > p2.Y && p.Y - t < p1.Y
        //    : p.Y + t > p1.Y && p.Y - t < p2.Y;

        double k, c;

        if (p2.X == p1.X)
        {
            return (p.X == p1.X && p.Y >= min(p1.Y, p2.Y) && p.Y <= max(p1.Y, p2.Y));
        }

        k = p2.X - p1.X == 0 ? 0 : ((double)(p2.Y - p1.Y)) / (p2.X - p1.X);

        if (k == 0)
        {
            return (p.Y == p1.Y && p.X >= min(p1.X, p2.X) && p.X <= max(p1.X, p2.X));
        }

        c = p1.Y - k * p1.X;

        return p.Y == p.X * k + c;
    }
};


#endif //STRAIGHT_RAY_EDGE_H
