#ifndef H_Spline3Interpolator_H
#define H_Spline3Interpolator_H

#include "DataPoints.h"
#include "TF1.h"
#include "TStyle.h"
#include <cmath>

#include "EqSolver.h"

class Spline3Interpolator : public DataPoints {

 public:
  Spline3Interpolator(int N=0, double *x=NULL, double *y=NULL, double alpha = 0, double beta = 0, TF1* fF0=NULL);
  Spline3Interpolator(vector<double>, vector<double>, double alpha = 0, double beta = 0, TF1* fF0=NULL);
  ~Spline3Interpolator() {
    // delete F0; // since F0 = fF0, there can only be ONE delete (cf. main.cpp)
    delete FInterpolator;
    delete FDerivator; // PROJECTO
    delete K;
  }

  double Interpolate(double x);
  void SetCurvatureLines(); // extra function
  void Draw(); //draw everything (points and interpolation function)
	TF1* GetInterpolationFunction() {
		return FInterpolator;
	}

  void SetFunction(TF1*);
  void Print(string FILE=""); // print results (optional) -- DONE

  /************ PROJECTO ***********/

  double Deriv(double x); // derivada de interpolação (PROJECTO)
  double fDerivator (double *fx, double *par){
    return Deriv(fx[0]);
  }

  TF1* GetDerivativeFunction() {
    return FDerivator;
  }

  /*********************************/

 private:
   friend int Binary_Search (int, double *, double);

   double fInterpolator(double *fx, double *par) {
     return Interpolate(fx[0]);
   }

  TF1* FInterpolator; //interpolation function
  TF1* F0;  //eventual underlying function
  double* K;
  TF1* FDerivator; // PROJECTO
  
};

#endif
