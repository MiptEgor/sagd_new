#ifndef FILT 
#define FILT
#include "mesh.h"
#include "const_values.h"
#include <iostream>
#include <iomanip>
#include <cmath>
class filtration 
{
public:
	filtration(mesh &Mesh, const_values cv);

	void process();
private:


	double k_perm(double theta);
	double kappa(double s);
	double D_kappa(double s);
	double lambda(double s);
	void calc_lay(double t);
	void find_max();
	void initial();


	mesh &area;
	const_values cv;
	double h, tau;
	double lambda_max;
	double max, min;
};

#endif //FILT