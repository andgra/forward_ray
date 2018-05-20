//
// Created by ANDGRA on 20.05.2018.
//

#ifndef BEAM_H
#define BEAM_H

#include "vec2.h"
#include "point.h"

class beam//луч
{
public:
    pointD startPoint;//координаты начальной точки
    vec2d direction;//направление
    double time;//время в начале луча
    double value;//значение амплитуды в начале луча
    int figureIndex;//номер фигуры
    //List<PointD> path;
    beam() {}
    beam(pointD startPoint, vec2d direction, double time, double value, int figureIndex)//, bool dif = false)
    {
        this->startPoint = startPoint;
        this->direction = vec2d::normalize(direction);
        this->time = time;
        this->value = value;
        this->figureIndex = figureIndex;
        //path = new List<PointD>();
        //path.Add(startPoint);
    }

    beam(beam oldBeam, pointD startPoint, vec2d direction, double time, double value,
         int figureIndex)//, bool reflected)
    {
        this->figureIndex = figureIndex;

        //this->path = new List<PointD>(oldBeam.path);
        //this->path.Add(startPoint);

        this->startPoint = startPoint;
        this->direction = direction;
        this->time = time;
        this->value = value;
    }
};

#endif //BEAM_H
