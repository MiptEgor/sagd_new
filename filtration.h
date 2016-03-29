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
	void calc_lay(double t);
	void initial();


	mesh &area;
	const_values cv;
	double h, tau;
	double lambda_max;
};

#endif //FILT