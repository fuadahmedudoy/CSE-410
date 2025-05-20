#include<bits/stdc++.h>
#include "Triangle.h"

Triangle::Triangle(point &p1, point &p2, point &p3)
{
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    rgb[0] = rand() % 256;
    rgb[1] = rand() % 256;
    rgb[2] = rand() % 256;
}

Triangle::Triangle(const Triangle &t)
{
    points[0] = t.points[0];
    points[1] = t.points[1];
    points[2] = t.points[2];
    rgb[0] = t.rgb[0];
    rgb[1] = t.rgb[1];
    rgb[2] = t.rgb[2];
}