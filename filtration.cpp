#include "filtration.h"
#include <iostream>

filtration::filtration(mesh &Mesh, const_values cv)
	:area(Mesh)
	, cv(cv)
	, h(cv.h)
	, tau(cv.tau)
	{}

double filtration::k_perm(double theta)
{
	return theta * theta;
}


void filtration::process()
{
	double t;
	int n = 0;
	std::cout<<cv.time_const/cv.tau<<std::endl;
	while (t<cv.time_const)
	{
		calc_lay();
		if (n % 100 == 0)
		{
			area.print_lay(n);
			std::cout<<n<<std::endl;
		}
		n++;
		t+=cv.tau;

	}


}


void filtration::calc_lay()
{
	area.get_left(0).W_l = 0;
	area.get_left(0).W_g = 0;

	for (int i = 0; i < area.get_n(); ++i)
	{
		area.get_right(i).W_l = - k_perm(area(i).theta_l()) / area(i).eta_l * (area(i).theta_s() * (cv.rho_s - cv.rho_l) 
			+ area(i).theta_g() * (cv.rho_g - cv.rho_l)) * cv.gravity;
		area.get_right(i).W_g = - k_perm(area(i).theta_l()) / area(i).eta_l * (area(i).theta_s() * (cv.rho_s - cv.rho_l) 
			+ area(i).theta_g() * (cv.rho_g - cv.rho_l)) * cv.gravity;

	}

	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).psi_l = area(i).psi_l - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau;
		area(i).psi_g = area(i).psi_g - (area.get_right(i).W_g - area.get_left(i).W_g) / h * tau;
	}




}