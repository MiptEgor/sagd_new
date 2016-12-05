#ifndef INTEGR 
#define INTEGR
#include "mesh.h"
#include "const_values.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "filtration.h"

class integrator
{
	mesh &area;
	filtration &Filtr;
	const_values cv;
public:
	integrator(mesh &Mesh, const_values CV, filtration &Filtration);
	void process();

};


#endif //INTEGR