#ifndef H_ODESOLVER_H
#define H_ODESOLVER_H
#include "ODEpoint.h"
#include "TFormula.h"

class ODEsolver{
  public:
    ODEsolver(){;}
    ODEsolver(vector<TFormula> Form);
    ~ODEsolver(){;};
    vector<ODEpoint> Eulersolver(const ODEpoint & P0, double xmin, double xmax, double h_step);
    vector<ODEpoint> RK2solver(const ODEpoint & P0, double xmin,double xmax, double h_step);
    vector<ODEpoint> RK4solver(const ODEpoint & P0, double xmin, double xmax, double h_step);
    vector<ODEpoint> RK4_AdapStep(const ODEpoint & P0, double xmin, double xmax, double h_step){;};
    vector<ODEpoint> Heun(const ODEpoint & P0, double xmin, double xmax, double h_step);
    //...
    void SetODEfunc(vector<TFormula> Form);
    vector<TFormula> GetODEfunc() const;

    void UpdateParameter(int, int, double); // PROJECTO
    ODEpoint RK4_iterator(const ODEpoint &, double step); // PROJECTO

  private:
    ODEpoint Heun_iterator (const ODEpoint &, double step);
    ODEpoint EULER_iterator (const ODEpoint &, double step);
    ODEpoint RK2_iterator (const ODEpoint &, double step);
    ODEpoint RK4_AS_iterator(const ODEpoint &, double step, vector<vector<double> > & K);
    // ODEpoint RK4_iterator(const ODEpoint &, double step);
    vector<TFormula> F;
};

#endif
