#ifndef MATRIX_H
#define MATRIX_H

class matrix
{
public:
    int row,col;
    double **mat;
    matrix(int row,int col);
    matrix(const matrix&);
    matrix& operator=(const matrix&);
    ~matrix();
};



#endif