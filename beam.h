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
    pointF startPoint;//координаты начальной точки
    vec2f direction;//направление
    float time;//время в начале луча
    float value;//значение амплитуды в начале луча
    int figureIndex;//номер фигуры
    //List<pointF> path;
    beam() {}
    beam(const pointF& startPoint, const vec2f& direction, float time, float value, int figureIndex)//, bool dif = false)
    {
        this->startPoint = startPoint;
        this->direction = vec2f::normalize(direction);
        this->time = time;
        this->value = value;
        this->figureIndex = figureIndex;
        //path = new List<pointF>();
        //path.Add(startPoint);
    }

    beam(const beam& oldBeam, const pointF& startPoint, const vec2f& direction, float time, float value,
         int figureIndex)//, bool reflected)
    {
        this->figureIndex = figureIndex;

        //this->path = new List<pointF>(oldBeam.path);
        //this->path.Add(startPoint);

        this->startPoint = startPoint;
        this->direction = direction;
        this->time = time;
        this->value = value;
    }
};

#endif //BEAM_H
