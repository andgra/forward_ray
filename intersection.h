//
// Created by ANDGRA on 20.05.2018.
//

#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "vec2.h"
#include "point.h"
#include "data.h"

class intersection {
public:
    pointD intersectionPoint;//точка пересечения
    double angle;//угол падения на грань к нормали
    vec2d direction;
    int figureIndex;//номер фигуры, к которой относится граница, в которую попал луч
    int edgeIndex;//номер фигуры, к которой относится граница, в которую попал луч
    double distance;
    double existed;

    vec2d normalVector;

    intersection() { existed = false; }

    intersection(pointD point, vec2d direction, int fIndex, int eIndex, vec2d normal, double distance) {
        intersectionPoint = point;
        this->direction = direction;
        this->figureIndex = fIndex;
        this->edgeIndex = eIndex;
        this->normalVector = normal;
        this->angle = acos(vec2d::dot(direction, normalVector) / (direction.length())) /
                     degreeToRadians;// * normalVector.length()));//длина нормального вектора всегда 1
        this->distance = distance;

        if (angle > 90)
            angle = fabs(angle - 180);
        existed = true;
    }
};

#endif //INTERSECTION_H
