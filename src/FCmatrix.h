#ifndef H_FCMATRIX_H
#define H_FCMATRIX_H
#include "Vec.h"

class FCmatrix{
  public:
    //constructors
    FCmatrix();
    FCmatrix(double** fM, int fm, int fn); // matrix fm xfn given from pointer of pointers
    FCmatrix(double* fM, int fm, int fn); // matrix fm x fn given as single pointer (what length ?!)
    FCmatrix(vector<Vec>); // matrix fm x fn given as vector of Vec


    // operators
    virtual Vec& operator[] (int) = 0; // get a row by giving índex inside []

    // methods
    virtual int Get_nRows() const = 0; // number of rows of M
    virtual int Get_nCols() const = 0; // number of columns of M
    virtual Vec GetRow(int i) const = 0; // retrieve row i
    virtual Vec GetCol(int i) const = 0; // retrieve column i
    virtual double Determinant() const = 0;

    // in row-i, return column-index of max element (in absolute value)
    virtual int GetRowMax(int i = 0) const = 0;
    // in column-j, return row-index (>=j) for which relative amplitude of Mij on the row is highest.
    virtual int GetColMax(int j = 0) const = 0; //

    virtual void Print() const; // print e.g. row by row (GetRow to the rescue...)
    vector<Vec> GetM() const; // virtual is not necessary here since we will not overload the function... (?)
    string GetClassname() const;

  protected:
    vector<Vec> M; // the matrix is a vector of Vec objects...
    string classname; // give the class a name...
};

#endif
