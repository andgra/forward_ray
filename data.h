//
// Created by ANDGRA on 20.05.2018.
//

#ifndef DATA_H
#define DATA_H

#include "helper.h"
#include <iostream>
#include <math.h>
#include "figure.h"
#include <algorithm>
#include <string>

using std::count;
using std::string;
using std::vector;

string fileName;
string filePath;
string fileDir;

string getOutPath() {
    return fileDir + "/out/";
}

figure background;
vector<figure> figureCollection;
int width;
int height;
int maxTime;

double dX = 0.0005;
double dY = 0.0005;
double dt = 0.000000036764;//время дискретизации по времени (0,000000036764)
double EPS = 1e-9;//значение, близкое к нулю
double degreeToRadians = M_PI / 180;//константа для перевода из градусов в радианы и обратно


double GetAbsorption(double coef, double distance)//рассчитываем поглощение
{
    //return (double)pow(exp(1), -coef * distance * dX * 100);//коэффициент в см^(-1), а расстояние в м
    return 1;//временно отключаем поглощение
}

bool isPointDif(pointD p, int figureIndex) {
    for(edge e: figureCollection[figureIndex].edges)
    {
        if (((round(e.xStart - p.X) == 0 && round(e.yStart - p.Y) == 0) ||
             (round(e.xEnd - p.X) == 0 && round(e.yEnd - p.Y) == 0)) && (figureIndex != 0)) {
            return true;
        }
    }
    return false;
}

int InWhichSmallestFigureIsPoint(pointD point) {
    double area = numeric_limits<double>::max();
    double curArea;
    int figureIndex = -1;
    for(figure f: figureCollection)
    {
        curArea = f.GetArea();
        if (f.ContainsPoint(point) && curArea < area) {
            figureIndex = f.index;
            area = curArea;
        }
    }
    return figureIndex;
}

int GetNextFigure(pointD p, vec2d dir, int coef = 1)//при 1-ном коэф.смотрим на 0,5 мм
{
    dir = vec2d::normalize(dir);
    return InWhichSmallestFigureIsPoint(pointD(p.X + dir.X * coef, p.Y + dir.Y * coef));
}

int GetNextFigure(pointD p) {
    auto dir = vec2d(0, -1);
    return InWhichSmallestFigureIsPoint(pointD(p.X + dir.X, p.Y + dir.Y));
}

int minSpeedCollection() {
    int min = numeric_limits<int>::max();
    for(figure f: figureCollection)
    {
        if(f.speed < min) min = f.speed;
    }
    return min;
}

#endif //DATA_H
