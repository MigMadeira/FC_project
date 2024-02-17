#ifndef H_VEC_H
#define H_VEC_H
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

class Vec{
  private:
    int N; // number of elements
    double *entries; // pointer to array of doubles

  public:
    Vec (int n = 1, double d = 0); // constructor with num. el. and value
    Vec (int, double *); // constructor with num. el. and ptr array
    Vec (vector<double>); // constructor vector<double>
    ~Vec(); // destructor
    Vec (const Vec &); // copy constructor...
    void SetEntries (int, double *); // set entries
    void Print(); // print the vector content
    void swap (int, int); // swap
    int size() const; // return size of vector
    double dot (const Vec &); // scalar product with another Vec
    Vec & operator= (const Vec &); // operator=
    Vec operator+ (const Vec &); // operator+
    Vec & operator+= (const Vec &); // operator+=
    Vec operator- (const Vec &); // operator-
    Vec & operator-= (const Vec &); // operator-=
    Vec operator+ () const; // operator+ (unary)
    Vec operator- () const; // operator- (unary)
    double & operator[] (int); // operator[]
    double operator[] (int) const; // operator[] when the Vec is a const
    Vec operator* (const Vec &) const; // operator* (CONST)
    Vec operator* (const double & scalar) const; // operator* a scalar (Vec LEFT)

    friend Vec operator* (const double & scalar, const Vec &); // operator* a scalar (Vec RIGHT)
};

#endif
