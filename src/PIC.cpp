#include "PIC.h"

int PIC::nshow = 0;

TCanvas PIC::Save_c ("Save_c", "PIC - Graph", 1500, 500);
TCanvas PIC::Save_dist ("Save_dist", "PIC - Distribution", 1000, 500);
TApplication PIC::App ("PIC", NULL, NULL);
TCanvas PIC::App_c ("App_c", "PIC - Real Time Plotting");

/****************** PLOT SETTINGS ***************************/

void PIC::SetGraphs(){

  phase_space = TGraph (npart);
  phase_space.SetTitle("Phase Space");
  phase_space.GetXaxis()->SetTitle("X");
  phase_space.GetXaxis()->SetLimits(x_min, x_max);
  phase_space.GetYaxis()->SetTitle("V");

  phase_space.SetMarkerColor(1);
  phase_space.SetMarkerSize(0.2);
  phase_space.SetMarkerStyle(8);

  poisson = TGraph (ngrid);
  poisson.SetTitle("Poisson");
  poisson.GetXaxis()->SetTitle("X");
  poisson.GetXaxis()->SetLimits(x_min, x_max);
  poisson.GetYaxis()->SetTitle("POT");

  // poisson.SetLineColor(2);
  // poisson.SetLineWidth(2);
  poisson.SetMarkerColor(1);
  poisson.SetMarkerSize(0.2);
  poisson.SetMarkerStyle(8);

  density = TGraph (ngrid);
  density.SetTitle("Density");
  density.GetXaxis()->SetTitle("X");
  density.GetXaxis()->SetLimits(x_min, x_max);
  density.GetYaxis()->SetTitle("DENS");

  // density.SetLineColor(2);
  // density.SetLineWidth(2);
  density.SetMarkerColor(1);
  density.SetMarkerSize(0.2);
  density.SetMarkerStyle(8);

}

void PIC::SetAppCanvas(){

  static TGraph tmp (1); // Template ...

  App_c.SetWindowSize(1500, 500);
  App_c.Divide(3,1);

  App_c.cd(1);
  App_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  tmp.SetTitle("Phase Space");
  tmp.GetXaxis()->SetTitle("X");
  tmp.GetXaxis()->SetLimits(x_min, x_max);
  tmp.GetYaxis()->SetTitle("V");
  tmp.Draw("AP");
  App_c.Update();
  App_c.Pad()->Draw();

  App_c.cd(2);
  App_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  tmp.SetTitle("Poisson");
  tmp.GetXaxis()->SetTitle("X");
  tmp.GetXaxis()->SetLimits(x_min, x_max);
  tmp.GetYaxis()->SetTitle("POT");
  tmp.Draw("AP");
  App_c.Update();
  App_c.Pad()->Draw();

  App_c.cd(3);
  App_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  tmp.SetTitle("Density");
  tmp.GetXaxis()->SetTitle("X");
  tmp.GetXaxis()->SetLimits(x_min, x_max);
  tmp.GetYaxis()->SetTitle("DENS");
  tmp.Draw("AP");
  App_c.Update();
  App_c.Pad()->Draw();

  gSystem->ProcessEvents();

}

void PIC::SetSaveCanvas(){

  // Save_c.SetWindowSize(1500, 500);
  Save_c.Divide(3,1);

  Save_c.cd(1);
  Save_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  // phase_space.Draw("AP");
  Save_c.Update();
  Save_c.Pad()->Draw();

  Save_c.cd(2);
  Save_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  // poisson.Draw("AP");
  Save_c.Update();
  Save_c.Pad()->Draw();

  Save_c.cd(3);
  Save_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);

  // density.Draw("AP");
  Save_c.Update();
  Save_c.Pad()->Draw();

}

void PIC::SetSaveDist(){

  // Save_dist.SetWindowSize(1000, 500);
  Save_dist.Divide(2,1);

  Save_dist.cd(1);
  Save_dist.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);
  Save_dist.Update();
  Save_dist.Pad()->Draw();


  Save_dist.cd(2);
  Save_dist.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);
  Save_dist.Update();
  Save_dist.Pad()->Draw();

}

/*********************************************************************/
/***************************** PIC **********************************/

PIC::PIC (vector<double> velocity, int Npart, double xmin, double xmax, int Ngrid) {
  // Ngrid é o número de pontos da grelha (logo, a caixa encontra-se dividida em (Ngrid - 1) subintervalos !!)
  vel_b = velocity;
  x_min = xmin;
  x_max = xmax;

  npart = Npart;

  ngrid = Ngrid; // ngrid >= 3 por causa do TridiagonalSolver()
  hgrid = (x_max - x_min) / ((double) (ngrid - 1));

  xgrid = new double [ngrid];
  dens_grid = new double [ngrid];
  pot_grid = new double [ngrid];

  for (int i = 0; i < ngrid; i++){
    xgrid[i] = x_min + hgrid * i;
  }

  // Reserve x_vpart memory (vector<ODEpoint>)
  x_vpart.reserve(npart);

  // ODEsolver Equation
  vector<TFormula> ODE;
  ODE.reserve(2);
  ODE.push_back(TFormula("dx/dt", "x[1]"));
  ODE.push_back(TFormula("dv/dt", "[0]")); // [0] == camp_el(r_i)
  ODE[1].SetParameter(0,0);

  Equation.SetODEfunc(ODE);

  // Building Tri_mat (ngrid-2 x ngrid-2) (constante)
  vector<double> d;
  d.reserve(ngrid - 3); // NGRID >= 3 !!

  vector<double> md;
  md.reserve(ngrid - 2);

  for (int i = 0; i < ngrid - 2; i++){
    md.push_back(-2);
  } // main diagonal

  for (int i = 0; i < ngrid - 3; i++){
    d.push_back(1);
  } // lower and upper diagonal

  vector<Vec> vector_mat;
  if (d.size() != 0) vector_mat.reserve(3);

  if (d.size() != 0) vector_mat.push_back(Vec(d));
  vector_mat.push_back(Vec(md));
  if (d.size() != 0) vector_mat.push_back(Vec(d));

  Tri_mat = FCmatrixBanded(vector_mat);

  // Densidade media (constante);
  dens_n0 = npart / (x_max - x_min);

  // Comprimento da caixa (constante)
  L = x_max - x_min;

  SetGraphs();
  SetSaveDist();
  SetSaveCanvas();

} // constructor


PIC::~PIC(){
  if (nshow != 0)
    App_c.WaitPrimitive();

  delete[] xgrid;
  delete[] pot_grid;
  delete[] dens_grid;
} // destructor


void PIC::FdistV (vector<double> veloc, bool save_plot, bool show_plot) {
  if (nshow != 0)
    App_c.WaitPrimitive();

  TRandom3 rgen(0);
  double * aux = new double [2];

  for (int i = 0; i < npart; i++){
    aux[0] = x_min + (x_max - x_min) * rgen.Uniform(0,1);

    if (rgen.Uniform(0,1) < 0.5)
      aux[1] = sqrt(-2*log(rgen.Uniform(0,1)))*cos(2*M_PI*rgen.Uniform(0,1)) + veloc[0];
    else
      aux[1] = sqrt(-2*log(rgen.Uniform(0,1)))*cos(2*M_PI*rgen.Uniform(0,1)) - veloc[1];

    x_vpart.push_back(ODEpoint(0, aux, 2));
  }

  delete[] aux;

  if (show_plot || save_plot){

    // v_inf e v_sup (indicativos) ...
    double v_inf = -5 * veloc[1];
    double v_sup = 5 * veloc[0];

    TF1 * func_dist = new TF1 ("Velocity Distribution Function", "1.0/(2.0*sqrt(2.0*pi))*(exp(-(x-[0])*(x-[0])/2.0) + exp(-(x+[1])*(x+[1])/2.0))", v_inf, v_sup);
    func_dist->SetParameters(veloc[0], veloc[1]);
    func_dist->SetNpx(5000);

    TH1F * hist_vel = new TH1F ("Vel_Hist", "Normalized Velocity Distribution", 1000, v_inf, v_sup);

    for (int i = 0; i < npart; i++)
      hist_vel->Fill(x_vpart[i][1]);

    double norm = hist_vel->GetEntries() * (v_sup - v_inf) / (double) (ngrid - 1);
    // double norm = hist_vel->Integral("width");
    hist_vel->Scale(1.0/norm);
    // cout << hist_vel->Integral("width") << endl;

    if (show_plot){

      App_c.Clear();
      App_c.SetWindowSize(1000, 500);
      App_c.Divide(2,1);

      App_c.cd(1);
      App_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);
      hist_vel->GetXaxis()->SetRangeUser(-2*veloc[1], 2*veloc[0]);
      hist_vel->Draw("HIST");
      App_c.Update();
      App_c.Pad()->Draw();

      App_c.cd(2);
      App_c.Pad()->SetMargin(0.12, 0.05, 0.1, 0.1);
      func_dist->SetTitle("Analytic Velocity Distribution");
      func_dist->GetXaxis()->SetRangeUser(-2*veloc[1], 2*veloc[0]);
      func_dist->Draw("C");
      App_c.Update();
      App_c.Pad()->Draw();

      App_c.WaitPrimitive();
      App_c.Clear();

      SetAppCanvas();
      nshow++;

    }

    if (save_plot){

      Save_dist.cd(1);
      hist_vel->GetXaxis()->SetRangeUser(-2*veloc[1], 2*veloc[0]);
      hist_vel->Draw("HIST");
      Save_dist.Update();
      Save_dist.Pad()->Draw();

      Save_dist.cd(2);
      func_dist->SetTitle("Analytic Velocity Distribution");
      func_dist->SetLineWidth(1);
      func_dist->GetXaxis()->SetRangeUser(-2*veloc[1], 2*veloc[0]);
      func_dist->Draw("C");
      Save_dist.Update();
      Save_dist.Pad()->Draw();

      Save_dist.Update();
      Save_dist.SaveAs("../Velocity_Distribution.eps");

    }

    delete hist_vel;
    delete func_dist;

  }

}

/****************** DENSIDADE DE PARTICULAS ************************/

int BinarySearch (int n, double * arr, double val){
  int left = 0;
  int right = n - 1;
  int middle = left + (right - left) / 2;

  /*
  if (val < arr[left] || val > arr[right]){
    cout << __PRETTY_FUNCTION__ << " : Search Err." << endl;
  }
  */

  while (left < right) {
    if (arr[middle] == val){
      return middle;
    }

    if (arr[middle] < val)
      left = middle + 1;

    else
      right = middle;

    middle = left + (right - left) / 2;

    // cout << left << "      " << right << "      " << middle << endl;
  }

  return left - 1;

} // return the LOWER closest index to 'val'. If 'val == x_max', the function returns (ngrid - 2) !!

void PIC::UpdateDensGrid(){
  // Para cada posição x da partícula, determinar a posição relativa na grelha e atualizar
  // os valores de densidade de partículas !!

  int j = 0;

  for (int i = 0; i < ngrid; i++){
    dens_grid[i] = 0;
  } // reset da densidade de partículas

  for (int i = 0; i < npart; i++){
    j = BinarySearch(ngrid, xgrid, x_vpart[i][0]);

    dens_grid[j] += (xgrid[j + 1] - x_vpart[i][0]) / (hgrid * hgrid);
    dens_grid[j + 1] += (x_vpart[i][0] - xgrid[j]) / (hgrid * hgrid);

    // if (xgrid[j] > x_vpart[i][0] || x_vpart[i][0] > xgrid[j + 1]) { cout << "NO - " << i << endl; break; }
  }

  // Condição Periódica (dens_grid[0] == dens_grid[ngrid - 1])
  double cst = dens_grid[0];
  dens_grid[0] += dens_grid[ngrid - 1];
  dens_grid[ngrid - 1] += cst;

  // Normalização
  for (int i = 0; i < ngrid; i++){
    dens_grid[i] = dens_grid[i] / dens_n0 - 1;
  }

/* METODO LINEAR (teste)

  double * test = new double [ngrid];
  for (int i = 0; i < ngrid; i++)
    test[i] = 0;

  auto start = high_resolution_clock::now();
  for (int i = 0; i < npart; i++){
    for (int k = 0; k < ngrid - 1; k++){
      if (xgrid[k] < x_vpart[i][0] && x_vpart[i][0] < xgrid[k + 1]){
        test[k] += (xgrid[k + 1] - x_vpart[i][0]) / (hgrid * hgrid);
        test[k + 1] += (x_vpart[i][0] - xgrid[k]) / (hgrid * hgrid);
        break;
      }
    }
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);
  cout << "Duration : " << duration.count() << endl;

  for (int i = 0; i < ngrid; i++){
    test[i] = test[i] / dens_n0 - 1;
    if (dens_grid[i] != test[i])
      cout << __PRETTY_FUNCTION__ << " : ERR - " << i << endl;;
  }

  delete[] test;
*/

} /* calcular a densidade de partículas ao longo da caixa (discretizada) */


void PIC::Density (bool save_plot, bool show_plot){

  if (nshow != 0)
    gSystem->ProcessEvents(); // Para interagir com a App ...

  UpdateDensGrid();

  if (show_plot || save_plot){
    for (int i = 0; i < ngrid; i++){
      density.SetPoint(i, xgrid[i], dens_grid[i]);
    }
  }

  if (show_plot){
    if (nshow == 0){
      SetAppCanvas();
      nshow++;
    }

    App_c.cd(3);
    density.GetXaxis()->SetLimits(x_min, x_max);
    density.Draw("AP");
    App_c.Update();
    App_c.Pad()->Draw();

    gSystem->ProcessEvents();

  }

  if (save_plot){

    Save_c.cd(3);
    density.GetXaxis()->SetLimits(x_min, x_max);
    density.Draw("AP");
    Save_c.Update();

    string Density_FileName = "../Density_" + to_string(x_vpart[0].Get_Time()).substr(0,5) + ".eps";
    Save_c.Pad()->SaveAs(Density_FileName.c_str());

  }

}

/****************************************************************/
/**************** POTENCIAL DAS PARTICULAS *********************/

void PIC::UpdatePotGrid(){

  Vec tri_b (ngrid - 2);
  for (int i = 0; i < ngrid - 2; i++)
    tri_b[i] = dens_grid[i + 1]*hgrid*hgrid;

  EqSolver Banded(Tri_mat, tri_b);

  // Solve the system ...
  Vec result;
  result = Banded.TridiagonalSolver();

  // Assign the private member pot_grid[] array pointer ...
  pot_grid[0] = 0;
  pot_grid[ngrid - 1] = 0;
  for (int i = 0; i < ngrid - 2; i++){
    pot_grid[i + 1] = result[i];
  }

} /* calcular o potencial discreto em cada ponto da grelha (i.e. em cada ponto em que
  a densidade também foi calculada) */

void PIC::Poisson (bool save_plot, bool show_plot){

  if (nshow != 0)
    gSystem->ProcessEvents(); // Para interagir com a App ...

  UpdatePotGrid();

  if (show_plot || save_plot){
    for (int i = 0; i < ngrid; i++){
      poisson.SetPoint(i, xgrid[i], pot_grid[i]);
    }
  }

  if (show_plot){
    if (nshow == 0){
      SetAppCanvas();
      nshow++;
    }

    App_c.cd(2);
    poisson.GetXaxis()->SetLimits(x_min, x_max);
    poisson.Draw("AP");
    App_c.Update();
    App_c.Pad()->Draw();

    gSystem->ProcessEvents();

  }

  if (save_plot){

    Save_c.cd(2);
    poisson.GetXaxis()->SetLimits(x_min, x_max);
    poisson.Draw("AP");
    Save_c.Update();

    string Poisson_FileName = "../Poisson_" + to_string(x_vpart[0].Get_Time()).substr(0,5) + ".eps";
    Save_c.Pad()->SaveAs(Poisson_FileName.c_str());
  }

}

/****************************************************************/
/*********************** TIME STEP *****************************/

void PIC::TimeStep (double dt){

  if (nshow != 0)
    gSystem->ProcessEvents(); // Para interagir com a App ...

  Spline3Interpolator Pot (ngrid, xgrid, pot_grid, dens_grid[0], dens_grid[ngrid - 1]); // dens_grid[0] == dens_grid[ngrid - 1]

  // cout << Pot.Deriv(x_min) << " -- E(x_min) ~ E(x_max) -- " << Pot.Deriv(x_max) << endl;

  for (int i = 0; i < npart; i++){
    Equation.UpdateParameter(1, 0, Pot.Deriv(x_vpart[i][0]));

    x_vpart[i] = Equation.RK4_iterator(x_vpart[i], dt);
    // x_vpart[i] = Equation.RK4solver(x_vpart[i], 0.0, dt, dt)[1]; -- Opção mais demorada !

    if (x_vpart[i][0] > x_max) x_vpart[i][0] = x_vpart[i][0] - L;
    if (x_vpart[i][0] < x_min) x_vpart[i][0] = x_vpart[i][0] + L;

  }
} /* avanço no tempo (por RK4_iterator) */


/*******************************************************************/
/******************** PHASE SPACE *********************************/

void PIC::Plot_Phase_Space (bool save_plot, bool show_plot){

  if (nshow != 0)
    gSystem->ProcessEvents(); // Para interagir com a App ...

  if (show_plot || save_plot){
    for (int i = 0; i < npart; i++){
      phase_space.SetPoint(i, x_vpart[i][0], x_vpart[i][1]);
    }
  }

  if (show_plot){
    if (nshow == 0){
      SetAppCanvas();
      nshow++;
    }

    App_c.cd(1);
    phase_space.GetXaxis()->SetLimits(x_min, x_max);
    phase_space.Draw("AP");
    App_c.Update();
    App_c.Pad()->Draw();

    gSystem->ProcessEvents();

  }

  if (save_plot){

    Save_c.cd(1);
    phase_space.GetXaxis()->SetLimits(x_min, x_max);
    phase_space.Draw("AP");
    Save_c.Update();

    string PhaseSpace_FileName = "../PhaseSpace_" + to_string(x_vpart[0].Get_Time()).substr(0,5) + ".eps";
    Save_c.Pad()->SaveAs(PhaseSpace_FileName.c_str());

  }
} /* Draw PhaseSpace */
