#include "const_values.h"
#include "mesh.h"
#include "filtration.h"
#include <stdio.h>
#include <iostream>


int main(int argc, char const *argv[])
{
	const_values CV;
	mesh Mesh(CV);
	Mesh(1).T = 273;
	filtration new_filtr(Mesh, CV);
	new_filtr.process();
	//std::cout <<Mesh(1).T<<std::endl;
	return 0;
}