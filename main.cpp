#include "const_values.h"
#include "mesh.h"
#include "filtration.h"
#include "integrator.h"
#include <stdio.h>
#include <iostream>


/*void process(const values)
	{
		double t = 0;
		int n = 0;
		std::cout<<std::setprecision(5);
		std::cout<<cv.time_const/cv.tau<<std::endl;
		while (t<cv.time_const)
		{
			if (n % 200 == 0)
			{
				area.print_lay(n);
				std::cout<<t<<std::endl;
			}
			calc_lay(t);
			n++;
			t+=tau;
			tau = h / 2 / lambda_max / 4;
			//if (tau * cv.W_bound_left / h > 0.5) tau = h / cv.W_bound_left / 2 / 2;
		}
	}*/

int main(int argc, char const *argv[])
{
	const_values CV;
	mesh Mesh(CV);
	Mesh(1).T = 273;
	filtration filtr_block(Mesh, CV);
	integrator integration_block(Mesh, CV, filtr_block);
	integration_block.process();
	//std::cout <<Mesh(1).T<<std::endl;
	return 0;
}