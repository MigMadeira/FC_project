#ifndef H_DataPoints_H
#define H_DataPoints_H

#include "TGraph.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <cstdio>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;

class DataPoints {
 public:
  DataPoints();
  DataPoints(int, double*, double*);
  DataPoints (vector<double>, vector<double>);
  virtual ~DataPoints();

  virtual double Interpolate(double x) {return 0.;}
  virtual void Draw();
  virtual void Print(string FILE="");
 protected:
  int N; // number of data points
  double *x, *y; // arrays
  static int Nplots;
};

#endif
