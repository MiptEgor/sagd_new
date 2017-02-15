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
	double calc_tau();
	void calc_lay(double t, double tau);
private:

	double viscosity(double T);
	double k_perm(double theta);
	
	void find_max();
	void initial();
	double W_ll(double l, double g);
	double W_lg(double l, double g);
	double W_gl(double l, double g);
	double W_gg(double l, double g);
	double first_lambda(double l, double g);
	double second_lambda(double l, double g);
	double W_l(double l, double g, double eta);
	double W_g(double l, double g, double eta);
	


	mesh &area;
	const_values cv;
	double h;
	double lambda_max;
	double max, min;
};

#endif //FILT