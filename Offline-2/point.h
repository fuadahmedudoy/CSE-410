#ifndef POINT_H
#define POINT_H
class point
{

    
public:
    double x,y,z,w;
    point();
    point(double x, double y, double z, double w);
    point(const point& other);
    ~point();
    point operator+(const point& p);
    point operator*(double scalar);
    friend std::ostream& operator<<(std::ostream&, const point&);
    
    
};
point cross(point p,point q);
double dot(point p,point q);
point rodrigues(point &, point &, double);

#endif