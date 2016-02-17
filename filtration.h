#ifndef FILT 
#define FILT
#include "mesh.h"
#include "const_values.h"
class filtration 
{
public:
	filtration(mesh &Mesh, const_values cv);

	void process();
private:


	double k_perm(double theta);
	void calc_lay(double t);
	void initial();

	const_values cv;
	mesh &area;
	double h, tau;
};

#endif //FILT