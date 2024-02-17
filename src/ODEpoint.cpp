#include "ODEpoint.h"

ODEpoint::ODEpoint(double tval, double* funct, int Ndimf){
  Ndim = Ndimf;
  t = tval;
  var.reserve(Ndimf);

  for (int i = 0; i < Ndimf; i++){
    var.push_back(funct[i]);
  }
} // constructor using double*

ODEpoint::ODEpoint(double tval, vector<double> funct){
  Ndim = funct.size();
  t = tval;
  var = funct;
} // constructor using vector<double>

ODEpoint::ODEpoint(const ODEpoint & copy){
  Ndim = copy.GetNdim();
  var = copy.Get_Var_vec();
  t = copy.Get_Time();
} // copy constructor

double ODEpoint::Get_Time() const{
  return t;
} // return only the independent variable (t)

vector<double> ODEpoint::Get_Var_vec() const{
  return var;
} // return the y1, ... , yNdim dependent variables

double * ODEpoint::Get_Var_ptr() const{
  double * Var_ptr = new double [Ndim];
  for (int i = 0; i < Ndim; i++)
    Var_ptr[i]=var[i];

  return Var_ptr;
} // same but as double *

double * ODEpoint::Get_VarTime() const{
  double * VarTime_ptr = new double [Ndim + 1];

  for (int i = 0; i < Ndim; i++)
    VarTime_ptr[i] = var[i];

  VarTime_ptr[Ndim] = t;
  return VarTime_ptr;
} // first the y1, ... , yNdim then t

int ODEpoint::GetNdim() const{
  return Ndim;
} // return the number of dependent variables

void ODEpoint::Set_Time(double tval){
  t = tval;
} // set independent variable

void ODEpoint::Set_Var(vector<double> funct){
  var = funct;
  Ndim = var.size();
} // set dependent variables

const double & ODEpoint::operator[](int i) const{
  return var[i];
}

double & ODEpoint::operator[] (int i){
  return var[i];
}

void ODEpoint::Print() const{
  cout << "y(t = " << t << ") = ";
  if (Ndim == 1){
    cout << var[0] << endl;
  }

  else{
    cout << "(";
    for (int i = 0; i < Ndim - 1; i++){
      cout << var[i] << " , " << flush;
    }
    cout << var[Ndim - 1] << ")" << endl;
  }
}

ODEpoint & ODEpoint::operator= (const ODEpoint & P){
  if (this != &P) {
    t = P.Get_Time();
    var = P.Get_Var_vec();
    Ndim = P.GetNdim();
  }

  return *this;
} // copy assignment
