#ifndef H_FCMATRIXBANDED_H
#define H_FCMATRIXBANDED_H
#include "FCmatrix.h"
#include "FCmatrixFull.h"

class FCmatrixBanded : public FCmatrix {
  public:
    // constructors
    FCmatrixBanded();
    FCmatrixBanded(double** fM, int fm); // matrix fM with fm columns
    FCmatrixBanded(double* fM, int fm); // "continuous" matrix fM with fm columns
    FCmatrixBanded(vector<Vec>); // vector<Vec> properly given as argument ...

    // copy constructor
    FCmatrixBanded (const FCmatrixBanded &);

    // operators
    FCmatrixBanded operator=(const FCmatrixBanded &); // equal 2 matrices of any kind (FCmatrix problem ...)
    FCmatrixBanded operator+(const FCmatrixBanded &) const; // add 2 matrices of any kind
    FCmatrixBanded operator-(const FCmatrixBanded &) const; // sub 2 matrices of any kind
    // Produtos de Tridiagonais sao Pentadiagonais (5 diagonais)
    // FCmatrixBanded operator*(const FCmatrixBanded &) const; // mul 2 matrices of any kind
    FCmatrixBanded operator*(double lambda) const; // mul matrix of any kind by scalar
    Vec operator*(const Vec &) const; // mul matrix by Vec

    // virtual inherited
    int Get_nRows() const; //number of rows of M
    int Get_nCols() const; // number of columns of M
    Vec GetRow(int i) const; // retrieve row i
    Vec GetCol(int i) const; // retrieve column i
    double Determinant() const{};
    void Print() const;

    Vec& operator[] (int);
    int GetRowMax(int i = 0) const{};
    int GetColMax(int j = 0) const{};
    // ...

    //extras
    Vec GetDiagonal(int i) const;
    friend FCmatrixFull convert (const FCmatrixBanded origin);
    friend FCmatrixBanded operator* (double, const FCmatrixBanded &);
};

#endif
