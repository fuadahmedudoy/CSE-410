#include<bits/stdc++.h>
#include"point.h"
using namespace std;

point::point(){
   this->x=0.0;
   this->y=0.0;
   this->z=0.0;
   this->w=1.0; 
}
point::point(double x, double y, double z, double w) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = w;
}
point::point(const point& other) {
    this->x = other.x;
    this->y = other.y;
    this->z = other.z;
    this->w = other.w;
}
point point::operator+(const point& q){
    double x=this->x+q.x;
    double y=this->y+q.y;
    double z=this->z+q.z;
    double w=1.0;
    return point(x, y, z, w);
}
point point::operator*(double scalar) {
    double x = this->x * scalar;
    double y = this->y * scalar;
    double z = this->z * scalar;
    double w = 1.0;

    return point(x, y, z, w);
}

point cross(point p,point q){
    double x=p.y*q.z -p.z*q.y;
    double y=p.z*q.y - p.x*q.z;
    double z=p.x*q.y - p.y*q.x;
    double w=1.0;
    return point(x,y,z,w);
}
double dot(point p,point q){
    return p.x*q.x +p.y*q.y +p.z*q.z;
}
std::ostream& operator<<(std::ostream& os, const point& p) {
    os << std::fixed << std::setprecision(7) << p.x << " " << p.y << " " << p.z << ' ';
    return os;
}

point rodrigues(point &x,point &axis,double angle){
    angle=angle*M_PI/180.0;
    double sine=sin(angle),cosine=cos(angle);
    double c=1.0-cosine;
    point p;
    p = x*cosine 
        + (axis * ( c * dot(x,axis))) 
        + (cross(axis,x)*sine);
    return p;
}

point::~point() {}