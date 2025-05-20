#ifndef TRIANGLE_H
#define TRIANGLE_H
#include "point.h"

class Triangle
{
public:
    point points[3];
    int rgb[3];
    Triangle(point&, point&, point&);
    Triangle(const Triangle&);
};

#endif