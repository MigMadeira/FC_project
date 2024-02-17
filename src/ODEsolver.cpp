#include "ODEsolver.h"

ODEsolver::ODEsolver(vector<TFormula> Form){
  F = Form;
} // constructor

void ODEsolver::SetODEfunc(vector<TFormula> Form){
  F = Form;
} // set vector<TFormula> (F)

void ODEsolver::UpdateParameter(int i, int p, double value){
  F[i].SetParameter(p, value);
} // set parameter 'p' from TFormula 'i' to 'value'

vector<TFormula> ODEsolver::GetODEfunc() const{
  return F;
} // return ODE system



/****************** EULER METHOD ************************/

ODEpoint ODEsolver::EULER_iterator (const ODEpoint & start, double step){
  int ydim = start.GetNdim();
  double* aux = start.Get_VarTime();

  vector<double> next_var;
  double next_t = start.Get_Time() + step;

  for (int i = 0; i < ydim; i++){
    next_var.push_back(start[i] + step * F[i].EvalPar(aux));
  }

  delete[] aux;
  return ODEpoint(next_t, next_var);
} // Euler Method Iterator

vector<ODEpoint> ODEsolver::Eulersolver(const ODEpoint & P0, double xmin, double xmax, double h_step){
  vector<ODEpoint> solution;
  int n = (xmax - xmin) / h_step;

  solution.push_back(P0);
  for (int i = 0; i < n; i++){
    solution.push_back(EULER_iterator(solution[i], h_step));
  }

  return solution;
} // Euler Method Solver


/****************** RK2 METHOD ************************/

ODEpoint ODEsolver::RK2_iterator (const ODEpoint & start, double step){
  int ydim = start.GetNdim();
  double * tandY = start.Get_VarTime();
  double * Y = start.Get_Var_ptr();
  double t = start.Get_Time();

  double * K1 = new double [ydim];
  double * K2 = new double [ydim];

  // Calculate K1
  for (int i = 0; i < ydim; i++){
    K1[i] = F[i].EvalPar(tandY);
  }

  // Calculate K2
  double * next_tandY = new double [ydim + 1];
  next_tandY[ydim] = t + step;
  for (int i = 0; i < ydim; i++){
    next_tandY[i] = Y[i] + step*K1[i];
  }

  for (int i = 0; i < ydim; i++){
    K2[i] = F[i].EvalPar(next_tandY);
  }

  // Now let's advance the y-array...
  for (int i = 0; i < ydim; i++){
    Y[i] += step/2.0 * (K1[i] + K2[i]);
  }

  t += step;
  ODEpoint next(t, Y, ydim);

  delete[] tandY;
  delete[] Y;
  delete[] K1;
  delete[] K2;
  delete[] next_tandY;

  return next;
} // RK2 Method Iterator

vector<ODEpoint> ODEsolver::RK2solver(const ODEpoint & P0, double xmin,double xmax, double h_step){
  vector<ODEpoint> solution;
  int n = (xmax - xmin) / h_step;

  solution.push_back(P0);
  for (int i = 0; i < n; i++){
    solution.push_back(RK2_iterator(solution[i], h_step));
  }

  return solution;
} // RK2 Method Solver


/****************** HEUN METHOD ************************/

ODEpoint ODEsolver::Heun_iterator (const ODEpoint & start, double step){
  int ydim = start.GetNdim();
  double * tandY = start.Get_VarTime();
  double * Y = start.Get_Var_ptr();
  double t = start.Get_Time();

  double * next_tandY = new double [ydim + 1];
  next_tandY[ydim] = t + step;
  for (int i = 0; i < ydim; i++)
    next_tandY[i] = Y[i] + step*F[i].EvalPar(tandY);

  for (int i = 0; i < ydim; i++)
    Y[i] += step/2.0 * (F[i].EvalPar(tandY) + F[i].EvalPar(next_tandY));

  t += step;
  ODEpoint next(t, Y, ydim);

  delete[] tandY;
  delete[] Y;
  delete[] next_tandY;

  return next;
} // Heun Method Iterator

vector<ODEpoint> ODEsolver::Heun(const ODEpoint & P0, double xmin, double xmax, double h_step){
  vector<ODEpoint> solution;
  int n = (xmax - xmin) / h_step;

  solution.push_back(P0);
  for (int i = 0; i < n; i++){
    solution.push_back(Heun_iterator(solution[i], h_step));
  }

  return solution;
} // Heun Method Solver


/********************* RK4 METHOD ***********************/

ODEpoint ODEsolver::RK4_iterator(const ODEpoint & start, double step){
  int ydim = start.GetNdim();
  double t = start.Get_Time();
  double * tandY = start.Get_VarTime();
  double * Y = start.Get_Var_ptr();

  double * K1 = new double [ydim];
  double * K2 = new double [ydim];
  double * K3 = new double [ydim];
  double * K4 = new double [ydim];

  double * next_tandY = new double [ydim + 1];

  // Calculate K1
  for (int i = 0; i < ydim; i++){
    K1[i] = F[i].EvalPar(tandY);
  }

  // Calculate K2
  next_tandY[ydim] = t + step / 2.0;
  for (int i = 0; i < ydim; i++){
    next_tandY[i] = Y[i] + step / 2.0 * K1[i];
  }

  for (int i = 0; i < ydim; i++){
    K2[i] = F[i].EvalPar(next_tandY);
  }

  // Calculate K3
  // next_tandY[ydim] = t + step / 2.0;
  for (int i = 0; i < ydim; i++){
    next_tandY[i] = Y[i] + step / 2.0 * K2[i];
  }

  for (int i = 0; i < ydim; i++){
    K3[i] = F[i].EvalPar(next_tandY);
  }

  // Calculate K4
  next_tandY[ydim] = t + step;
  for (int i = 0; i < ydim; i++){
    next_tandY[i] = Y[i] + step * K3[i];
  }

  for (int i = 0; i < ydim; i++){
    K4[i] = F[i].EvalPar(next_tandY);
  }

  // Calculate next Y and t
  for (int i = 0; i < ydim; i++)
    Y[i] += step/6.0 * (K1[i] + K2[i] + K2[i] + K3[i] + K3[i] + K4[i]);

  t += step;
  ODEpoint next(t, Y, ydim);

  delete[] tandY;
  delete[] Y;
  delete[] next_tandY;
  delete[] K1;
  delete[] K2;
  delete[] K3;
  delete[] K4;

  return next;

} // RK4 Method Iterator

vector<ODEpoint> ODEsolver::RK4solver(const ODEpoint & P0, double xmin, double xmax, double h_step){
  vector<ODEpoint> solution;
  int n = (xmax - xmin) / h_step;
  solution.reserve(n + 1);

  solution.push_back(P0);
  for (int i = 0; i < n; i++){
    solution.push_back(RK4_iterator(solution[i], h_step));
  }

  return solution;
} // RK4 Method Solver
