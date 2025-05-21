#include<bits/stdc++.h>
#include "Triangle.h"


static unsigned long int g_seed = 1;
inline int genrand()
{
g_seed = (214013 * g_seed + 2531011);
return (g_seed >> 16) & 0x7FFF;
}
Triangle::Triangle(point &p1, point &p2, point &p3)
{
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    rgb[0] = genrand() % 256;
    rgb[1] = genrand() % 256;
    rgb[2] = genrand() % 256;
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