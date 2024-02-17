#include "Spline3Interpolator.h"

using namespace std;

int Binary_Search (int n, double * arr, double val){
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

    // cout << left << "     " << right << "      " << middle << endl;
  }

  return left - 1;
} // return the LOWER closest index to 'val'

Spline3Interpolator::Spline3Interpolator(int fN, double *fx, double *fy, double alpha, double beta, TF1* fF0) : DataPoints(fN, fx, fy) {
  // DataPoints::Print();
  F0=fF0;

	/*  In order to define a function coming from somewhere in my Classes i need to
		use a special constructor for TF1 (see https://root.cern.ch/doc/master/exampleFunctor_8C.html) i.e.

	  TF1 * f3 = new TF1("f3",intg,&MyIntegFunc::Integral, xmin, xmax, 0);

		only in our current case the function is defined in the class as Spline3Interpolator::fInterpolator
		and the object that holds such a fInterpolator method is precisely the object i am at (((-;
		thus the "this" applies here perfectly...
			Now, to understand that it DOES work, note that fInterpolator actually does this

			return Interpolate(fx[0]);

		meaning it evaluates the Interpolate method with a generic input. It is precisely the
		Interpolate method that actually implements/evaluates the interpolation algorithm....

	*/

  FInterpolator = new TF1("FInterpolator", this, &Spline3Interpolator::fInterpolator, x[0]-0.001 ,x[N-1]+0.001, 0);
  FDerivator = new TF1 ("FInterpolator", this, &Spline3Interpolator::fDerivator, x[0]-0.001 , x[N - 1]+0.001, 0); // PROJECTO

  K = new double [N];
  K[0] = alpha;
  K[N - 1] = beta;
  SetCurvatureLines();
}

Spline3Interpolator::Spline3Interpolator(vector<double> fx, vector<double> fy, double alpha, double beta, TF1* fF0) : DataPoints(fx, fy) {
  // DataPoints::Print();
  F0=fF0;

  FInterpolator = new TF1("FInterpolator", this, &Spline3Interpolator::fInterpolator, x[0]-0.001 ,x[N-1]+0.001, 0);
  K = new double [N];
  K[0] = alpha;
  K[N - 1] = beta;
  SetCurvatureLines();
}


void Spline3Interpolator::SetCurvatureLines() {
  // Building Tri_mat N-2 x N-2
  vector<double> d;
  d.reserve(N - 3);

  vector<double> md;
  md.reserve(N - 2);

  Vec tri_b(N - 2);

  for (int i = 0; i < N - 2; i++){
    md.push_back(2*(x[i] - x[i + 2]));

    // Building tri_b
    tri_b[i] = 6 * (((y[i] - y[i + 1])/(x[i] - x[i + 1])) - ((y[i + 1] - y[i + 2])/(x[i + 1] - x[i + 2])));
  } // main diagonal

  tri_b[0] = tri_b[0] - K[0]*(x[0] - x[1]);
  tri_b[N - 3] = tri_b[N - 3] - K[N - 1] * (x[N - 2] - x[N - 1]);

  for (int i = 1; i < N - 2; i++){
    d.push_back(x[i] - x[i + 1]);
  } // lower and upper diagonal

  vector<Vec> vector_mat;

  if (d.size() != 0) vector_mat.push_back(Vec(d));
  vector_mat.push_back(Vec(md));
  if (d.size() != 0) vector_mat.push_back(Vec(d));

  FCmatrixBanded Tri_mat(vector_mat);
  EqSolver banded(Tri_mat,tri_b);
  // Solve the system....
  Vec result;
  result = banded.TridiagonalSolver();

  // Assign the private member K[] array pointer ...
  for (int i = 0; i < N - 2; i++){
    K[i + 1] = result[i];
  }
}

double Spline3Interpolator::Interpolate(double fx) {
  double A;
  int i = 0;

  // detect in which segment is x, if outside flag as extrapolation !
  if (x[0] <= fx && fx <= x[N - 1])
    i = Binary_Search(N, x, fx);

  if (fx < x[0]){
    cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " < x[0] = " << x[0 ] << ") " << endl;
    i = 0;
  }

  if (fx > x[N - 1]){
    cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " > x[N - 1] = " << x[N - 1] << ") " << endl;
    i = N - 2;
  }

  /*
  // detect in which segment is x, if outside flag as extrapolation !
  for (int j = 0; j < N - 1; j++){
    if (fx < x[0]){
      cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " < x[0] = " << x[0 ] << ") " << endl;
      i = 0;
      break;
    }

    if (fx > x[N - 1]){
      cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " > x[N - 1] = " << x[N - 1] << ") " << endl;
      i = N - 2;
      break;
    }

    if (x[j] <= fx && fx < x[j + 1]){
      i = j;
      break;
    }
  }
  */

  // it is pointless to use all the spline functions... we need just one !
  //...

  A = 1./6. * K[i] * ((pow(fx - x[i + 1], 3))/(x[i] - x[i + 1]) - (fx - x[i + 1])*(x[i] - x[i + 1])) -
      1./6. * K[i + 1] * ((pow(fx - x[i], 3))/(x[i] - x[i + 1]) - (fx - x[i])*(x[i] - x[i + 1])) +
      (y[i]*(fx - x[i + 1])  - y[i + 1]*(fx - x[i]))/(x[i] - x[i + 1]);

  // cout<<endl;
  // cout<<"Spline3["<<xval<<"]="<<A<<flush;

  return A;
}

void Spline3Interpolator::Draw() {
  if (Nplots == 0) {
    //create application
    TApplication * MyRootApp;
    MyRootApp = new TApplication("Spline3Interpolator", NULL, NULL);
    MyRootApp->SetReturnFromRun(true);
  }
  TCanvas *c2 = new TCanvas("c1","c1",600,500);

  gPad->Clear();
  gStyle->SetLabelSize(0.06,"X");
  gStyle->SetLabelSize(0.06,"Y");
  TGraph *g = new TGraph(N,x,y);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kGreen);
  g->SetMarkerSize(2.5);

	FInterpolator->SetLineColor(kRed);
	FInterpolator->SetLineWidth(4);
	FInterpolator->SetTitle("Spline3");
  if (F0 != NULL) {
    F0->SetLineColor(kBlue);
    F0->SetLineWidth(4);
    F0->SetLineStyle(2);
  }

  gPad->DrawFrame(FInterpolator->GetXmin(), FInterpolator->GetMinimum(),
                  FInterpolator->GetXmax(), FInterpolator->GetMaximum(), FInterpolator->GetTitle());
  g->Draw("PSAME");
  FInterpolator->SetNpx(10000); // Necessary for the Pendulum Experimental Data (LMOO) Interpolation (cf. main.cpp)
  FInterpolator->Draw("SAME");
  if (F0 != NULL)
    F0->Draw("SAME");

  c2->Update();
  c2->SaveAs("plot_Spline3.eps");
  gPad->WaitPrimitive();
  // delete MyRootApp;
  delete g;
  delete c2;
  Nplots++;
}

void Spline3Interpolator::SetFunction(TF1* f) {
  F0 = f;
}

void Spline3Interpolator::Print(string FILE){
  int Npoints = 5000; // can be changed if necessary ...
  double x_results;

  printf("[Spline3::Print]: ");
  for (int i = 0; i <= Npoints; i++){
    x_results = (double) i / (double) Npoints *(x[N - 1] - x[0]);
    printf("[Spline3] (x=%f, y=%f) ", x_results, FInterpolator->Eval(x_results));
  }

  printf("\n");
  if (FILE != "")
  {
    printf("FILE: %s\n", FILE.c_str());
    double* xr = new double [Npoints];
    double* yr = new double [Npoints];

    ofstream outfile (FILE + ".txt");

  	for (int n = 0; n < Npoints; n++) {
      xr[n] = (double) n / (double) Npoints *(x[N - 1] - x[0]);
      yr[n] = FInterpolator->Eval(xr[n]);
  		outfile << setprecision(9);
  		outfile << fixed;
  		outfile << setw(17) << "[Spline3] :" <<
  		setw(20) << xr[n] << setw(20) << yr[n] << endl;
  	}

  	outfile.close();

    /******* CREATE GRAPH WITH SPLINE3 AS RED LINE AND DATAPOINTS AS ... POINTS *************/

    TCanvas* cr = new TCanvas ("cr", "Spline3 Results", 500, 500);
    cr->SetGrid();

    TGraph* gr1 = new TGraph (Npoints, xr, yr);
    gr1->SetTitle(FILE.c_str());
    gr1->GetXaxis()->SetTitle("X Axis");
    gr1->GetXaxis()->SetTitleOffset(1.4);
    gr1->GetYaxis()->SetTitle("Y Axis");
    gr1->GetYaxis()->SetTitleOffset(1.2);
    /* gr1->SetMarkerColor(1);
    gr1->SetMarkerSize(0.1);
    gr1->SetMarkerStyle(4); */
    gr1->SetLineColor(2);
    gr1->SetLineWidth(1.1);
    gr1->Draw("AC");

    TGraph* gr2 = new TGraph (N, x, y);
    gr2->SetMarkerColor(1);
    gr2->SetMarkerSize(0.2);
    gr2->SetMarkerStyle(8);
    gr2->Draw("P");

    cr->Modified();
    cr->Update();
    cr->SaveAs((FILE+".eps").c_str());

    delete gr1, cr, gr2, xr, yr;
  }
};


/************* PROJECTO - DERIVADA (INTERPOLAÇÃO) *****************/

double Spline3Interpolator::Deriv(double fx){
  double Derivative;
  int i = 0;

  // detect in which segment is x, if outside flag as extrapolation !
  if (x[0] <= fx && fx <= x[N - 1]){
    i = Binary_Search(N, x, fx);
    if (x[i] > fx || fx > x[i + 1])
      cout << "NOOOO" << endl;
  }

  if (fx < x[0]){
    cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " < x[0] = " << x[0 ] << ") " << endl;
    i = 0;
  }

  if (fx > x[N - 1]){
    cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " > x[N - 1] = " << x[N - 1] << ") " << endl;
    i = N - 2;
  }

  /*
  // detect in which segment is x, if outside flag as extrapolation !
  for (int j = 0; j < N - 1; j++){
    if (fx < x[0]){
      cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " < x[0] = " << x[0 ] << ") " << endl;
      i = 0;
      break;
    }

    if (fx > x[N - 1]){
      cout << __PRETTY_FUNCTION__ << " : Extrapolation (" << fx << " > x[N - 1] = " << x[N - 1] << ") " << endl;
      i = N - 2;
      break;
    }

    if (x[j] <= fx && fx < x[j + 1]){
      i = j;
      break;
    }
  }
  */

  Derivative = 1./6. * K[i] * (((3.0*pow(fx - x[i + 1], 2)) / (x[i] - x[i + 1])) - (x[i] - x[i + 1])) -
               1./6. * K[i + 1] * (((3.0*pow(fx - x[i], 2)) / (x[i] - x[i + 1])) - (x[i] - x[i + 1])) +
               ((y[i] - y[i + 1]) / (x[i] - x[i + 1]));

  // cout<<endl;
  // cout<<"Spline3Derivative["<<fx<<"]="<<Derivative<<flush;
  return Derivative;
} // derivada em qualquer ponto do dominio de interpolaçao
