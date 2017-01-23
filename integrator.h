#ifndef INTEGR 
#define INTEGR
#include "mesh.h"
#include "const_values.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include "filtration.h"
#include "energy.h"

class integrator
{
	mesh &area;
	filtration &Filtr;
	energy &energy_block;
	const_values cv;
public:
	integrator(mesh &Mesh, const_values CV, filtration &Filtration, energy &Energy);
	void process();

};


#endif //INTEGR