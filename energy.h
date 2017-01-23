#ifndef ENERGY
#define ENERGY
#include "const_values.h"
#include "mesh.h"

class energy
{
public:
	energy(mesh &Mesh, const_values CV);
	void calc_lay(double t, double tau);

private:
	mesh &area;
	const_values cv;
	double h;
};






#endif //ENERGY