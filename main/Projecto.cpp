#include <iostream>
#include "PIC.h"
// ...

using namespace std;

int main () {

	bool save = false;
	bool show = false;
	bool show_real_time = true;
	int n = 500;

	double xmin = 0.0;
	double xmax = 50.0;
	int Npart = 100000;
	int Ngrid = 1001;
	double dt = 0.04;
	vector<double> velvec(2);
	velvec[0] = 5; //velocity of first beam (negative)
	velvec[1] = 5; //velocity of second beam (positive)

	PIC twofluid = PIC(velvec,Npart,xmin,xmax,Ngrid);

	twofluid.FdistV(velvec, save, false);
	twofluid.Density(save, show);
	twofluid.Poisson(save, show);
	twofluid.Plot_Phase_Space(save, show);

	auto start = high_resolution_clock::now();

	for (int i = 1; i <= 1500; i++) {
		if (i % n == 0 && i >= n) save = false;

		twofluid.TimeStep(dt);
		twofluid.Density(save, show_real_time);
		twofluid.Poisson(save, show_real_time);
		twofluid.Plot_Phase_Space(save, show_real_time);

		// cout << i << endl;

		if (i % n == 0 && i >= n) save = false;
	}

	auto stop = high_resolution_clock::now();
	auto duration = duration_cast<seconds>(stop - start);
	cout << "Duration : " << duration.count() << endl;

	return 0;
}
