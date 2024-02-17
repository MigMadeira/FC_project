#include "DataPoints.h"

int DataPoints::Nplots=0;

DataPoints::DataPoints() {
  N = 0;
  x = NULL;
  y = NULL;
}

DataPoints::DataPoints(int fN, double* fx, double* fy) : N(fN) {
  x = new double[N];
  y = new double[N];

  for (int i=0; i<N; i++) {
    x[i] = fx[i];
    y[i] = fy[i];
    // printf("[Datapoints] x=%f, y=%f \n", x[i], y[i]);
  }
}

DataPoints::DataPoints (vector<double> fx, vector<double> fy) {
  N = fx.size();
  x = new double [N];
  y = new double [N];

  for (int i=0; i<N; i++){
    x[i] = fx[i];
    y[i] = fy[i];
    // printf("[Datapoints] x=%f, y=%f \n", x[i], y[i]);
  }
}

DataPoints::~DataPoints() {
  if (x != NULL || y != NULL){
    delete [] x;
    delete [] y;
  }
}

void DataPoints::Draw() {
  TGraph *g = new TGraph(N,x,y);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kRed);
  g->SetMarkerSize(2.5);
  if (Nplots == 0) {
    //create application
    TApplication * MyRootApp;
    MyRootApp = new TApplication("DataPoints", NULL, NULL);
    MyRootApp->SetReturnFromRun(true);
  }
  TCanvas *c0 = new TCanvas("c0","c0",600,500);
  g->Draw("PA");
  c0->Update();
  gPad->WaitPrimitive();
  delete g;
  Nplots++;
}

void DataPoints::Print(string FILE) {

  printf("[Datapoints::Print]: ");
  for (int i=0; i<N; i++) {
    printf("[Datapoints] (x=%f, y=%f) ", x[i], y[i]);
  }
  printf("\n");
  if (FILE != "")
  {
    printf("FILE: %s\n", FILE.c_str());

    ofstream outfile (FILE + ".txt");

  	for (int n = 0; n < N; n++) {
  		outfile << setprecision(9);
  		outfile << fixed;
  		outfile << setw(17) << "[Datapoints] :" <<
  		setw(20) << x[n] << setw(20) << y[n] << endl;
  	}

  	outfile.close();
  }
}
