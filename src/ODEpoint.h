#ifndef H_ODEPOINT_H
#define H_ODEPOINT_H
#include <iostream>
#include <vector>
using namespace std;

class ODEpoint {
  public:
    ODEpoint() {t = 0; Ndim = 0;}; // default constructor
    ~ODEpoint(){;} // destructor
    ODEpoint(double, double*, int); // using double *
    ODEpoint(double, vector<double>); // using vector<double>
    ODEpoint(const ODEpoint &); // copy constructor

    // Member Access Functions
    vector<double> Get_Var_vec() const; // return the y1, ... , yNdim dependent variables
    double * Get_Var_ptr() const; // same but as double *
    double * Get_VarTime() const; // first the y1, ... , yNdim then t
    int GetNdim() const; // return the number of dependent variables
    double Get_Time() const; // return only the independent variable (t)

    void Set_Time(double); // set independent variable
    void Set_Var(vector<double>); // set dependent variables

    // Operators
    ODEpoint & operator= (const ODEpoint &); // copy assignment
    const double & operator[] (int) const;
    double & operator[] (int);
    void Print() const;

  private:
    double t;
    vector<double> var;
    int Ndim;
};

#endif
