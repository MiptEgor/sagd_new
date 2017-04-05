#include "const_values.h"
#include "mesh.h"
#include "filtration.h"
#include "energy.h"
#include "integrator.h"
#include <stdio.h>
#include <iostream>
 #include <fenv.h>



int main(int argc, char const *argv[])
{
	//feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);
	const_values CV;
	mesh Mesh(CV);
	filtration filtr_block(Mesh, CV);
	energy energy_block(Mesh, CV);
	integrator integration_block(Mesh, CV, filtr_block, energy_block);
	integration_block.process();
	//std::cout <<Mesh(1).T<<std::endl;
	return 0;
}