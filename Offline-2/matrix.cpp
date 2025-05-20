#include<bits/stdc++.h>
#include"matrix.h"
#include"point.h"

matrix::matrix(int row,int col){
    this->row=row;
    this->col=col;
    this->mat=new double*[row];
    for(int i=0;i<row;i++){
        this->mat[i]=new double[col];
        for(int j=0;j<col;j++){
            this->mat[i][j]=(i==j) ? 1.0 : 0.0;
        }
    }
}

matrix::matrix(const matrix& m){
    this->row=m.row;
    this->col=m.col;
    this->mat=new double*[row];
    for(int i=0;i<row;i++){
        this->mat[i]=new double[col];
        for(int j=0;j<col;j++){
            this->mat[i][j]=m.mat[i][j];
        }
    }
}

matrix::~matrix() {
    for (int i = 0; i < this->row; i++) {
        delete[] this->mat[i];
    }
    delete[] this->mat;
}

matrix& matrix::operator=(const matrix& other) {
    if (this == &other) {
        return *this;
    }

    for (int i = 0; i < this->row; i++) {
        delete[] this->mat[i];
    }
    delete[] this->mat;

    this->row = other.row;
    this->col = other.col;
    this->mat = new double*[this->row];

    for (int i = 0; i < this->row; i++) {
        this->mat[i] = new double[this->col];
        for (int j = 0; j < this->col; j++) {
            this->mat[i][j] = other.mat[i][j];
        }
    }

    return *this;
}
