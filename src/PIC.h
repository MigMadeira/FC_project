#ifndef H_PIC_H
#define H_PIC_H

#include <iostream>
#include <vector>
#include <chrono>
#include "TRandom3.h"
#include "ODEpoint.h"
#include "ODEsolver.h"
#include "EqSolver.h"
#include "Spline3Interpolator.h"
#include "TH1F.h"
#include "TApplication.h"
#include "TSystem.h"

using namespace std;
using namespace chrono;

class PIC {
  public :
    PIC () {;}; // default constructor
    PIC (vector<double> velocity, int Npart = 1000, double xmin = 0.0, double xmax = 1.0, int Ngrid = 100); // constructor
    ~PIC(); // destructor

    void FdistV (vector<double> veloc, bool save_plot, bool show_plot = false);
    void Plot_Phase_Space (bool save_plot, bool show_plot = false);
    void Density (bool save_plot, bool show_plot = false);
    void Poisson (bool save_plot, bool show_plot = false);
    void TimeStep (double dt);

  protected :
    static int nshow;

    static TCanvas Save_dist;
    static TCanvas Save_c;
    static TApplication App;
    static TCanvas App_c;

  private :
    vector<double> vel_b; // vector de 2 velocidades (vb), ambas positivas, a subtrair e somar na distribuição de velocidades
    double x_min; // equivalente ao xmin do enunciado
    double x_max; // equivalente ao xmax do enunciado

    int npart; // número de partículas usada na simulação

    int ngrid; // número de divisões da grelha
    double hgrid; // dx

    double * xgrid; // array com as 'ngrid' posições da grelha considerada (onde são calculados o potencial e a densidade de particulas)
    double * dens_grid; // array de 'ngrid' doubles com a densidade nas posições da grelha considerada
    double * pot_grid; // array de 'ngrid' doubles com o potencial nas posições da grelha considerada

    vector<ODEpoint> x_vpart; // vector com a posição e velocidade de cada partícula-i (num dado instante 't')
    ODEsolver Equation; // objeto da classe ODEsolver onde são resolvidas as equações do nosso sistema

    FCmatrixBanded Tri_mat; // matriz tridiagonal usada no calculo do potencial
    double dens_n0; // densidade média de partículas por unidade de comprimento
    double L; // comprimento da caixa

    TGraph phase_space; // gráfico do espaço de fases
    TGraph poisson; // gráfico do potencial
    TGraph density; // gráfico da densidade de partículas

    void SetAppCanvas();
    void SetGraphs();
    void SetSaveCanvas();
    void SetSaveDist();

    void UpdateDensGrid();
    void UpdatePotGrid();

    friend int BinarySearch (int, double *, double);

};

#endif
