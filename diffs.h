//
// Created by ANDGRA on 21.05.2018.
//

#ifndef DIFFS_H
#define DIFFS_H

#include "point.h"
#include <unordered_map>

typedef std::unordered_map<int, double> diffMap;
typedef std::unordered_map<pointD*, diffMap> diffsMap;

class diffs {
public:
    diffsMap inst;
    diffs() {}

};
#endif //DIFFS_H
