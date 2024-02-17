#ifndef H_FCMATRIXFULL_H
#define H_FCMATRIXFULL_H
#include "FCmatrix.h"
#include <cmath>

class FCmatrixFull : public FCmatrix {
  public:
    // constructors
    FCmatrixFull();
    FCmatrixFull(double** fM, int fm, int fn); //matrix fm x fn
    FCmatrixFull(double* fM, int fm, int fn);
    FCmatrixFull(vector<Vec>);

    // copy constructor
    FCmatrixFull (const FCmatrixFull&);

    // operators
    FCmatrixFull operator=(const FCmatrix &); // equal 2 matrices of any kind
    FCmatrixFull operator+(const FCmatrix &) const; // add 2 matrices of any kind
    FCmatrixFull operator-(const FCmatrix &) const; // sub 2 matrices of any kind
    FCmatrixFull operator*(const FCmatrix &) const; // mul 2 matrices of any kind
    FCmatrixFull operator*(double lambda) const; // mul matrix of any kind by scalar
    friend FCmatrixFull operator* (double, const FCmatrixFull &);

    Vec operator*(const Vec &) const; // mul matrix by Vec

    // virtual inherited
    int Get_nRows() const; //number of rows of M
    int Get_nCols() const; // number of columns of M
    Vec GetRow(int i) const; // retrieve row i
    Vec GetCol(int i) const; // retrieve column i
    double Determinant() const;

    Vec& operator[] (int);
    int GetRowMax(int i = 0) const;
    int GetColMax(int j = 0) const;
    // ...
    void swapRows(int, int);
};

#endif
