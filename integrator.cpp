#include "integrator.h"

integrator::integrator(mesh &Mesh, const_values CV, filtration &Filtration)
	: area(Mesh)
	, Filtr(Filtration)
	, cv(CV)
	{};

void integrator::process()
{
	double tau = cv.tau;
	double t = 0;
	int n = 0;
	std::cout<<std::setprecision(5);
	//std::cout<<cv.time_const/cv.tau<<std::endl;
	while (t<cv.time_const)
	{
		if (n % 200 == 0)
		{
			area.print_lay(n);
			std::cout<<t<<std::endl;
		}
		Filtr.calc_lay(t, tau);
		n++;
		t+=tau;

		tau = Filtr.calc_tau();
		//if (tau * cv.W_bound_left / h > 0.5) tau = h / cv.W_bound_left / 2 / 2;
	}
}