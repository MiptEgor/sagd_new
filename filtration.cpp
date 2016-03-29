#include "filtration.h"


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

void filtration::initial()
{
	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).psi_l = 1.;
		area(i).psi_g = 0.;
	}

}

void filtration::process()
{
	initial();
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
		tau = h / 2 / lambda_max;
	}
}


void filtration::calc_lay(double t)
{
	area.get_left(0).W_l = 0;
	area.get_left(0).W_g = 4e-6;
	double a = 0;
	lambda_max = 0;
	for (int i = 0; i < area.get_n(); ++i)
	{
		area.get_right(i).W_l =  cv.K_abs * k_perm(area(i).theta_l()) / area(i).eta_l * (area(i).theta_s() * (cv.rho_s - cv.rho_l) 
			+ area(i).theta_g() * (cv.rho_g - cv.rho_l)) * cv.gravity;
		area.get_right(i).W_g =  cv.K_abs * k_perm(area(i).theta_g()) / area(i).eta_g * (area(i).theta_s() * (cv.rho_s - cv.rho_g) 
			+ area(i).theta_l() * (cv.rho_l - cv.rho_g)) * cv.gravity;
	//	if (area.get_right(i).W_l > W_max) W_max = area.get_right(i).W_l;
	//	if (area.get_right(i).W_g > W_max) W_max = area.get_right(i).W_g;
	}

	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).psi_l = area(i).psi_l - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau;
		area(i).psi_g = area(i).psi_g - (area.get_right(i).W_g - area.get_left(i).W_g) / h * tau;

		//Определение максимального собственного значения
		a = fabs(2*cv.K_abs / cv.eta_l * (area(i).theta_g() * (cv.rho_g - cv.rho_l) + area(i).theta_s() * (cv.rho_s - cv.rho_l)));
		if (lambda_max < a) lambda_max = a;
		a = fabs(2*cv.K_abs / cv.eta_g * (area(i).theta_l() * (cv.rho_l - cv.rho_g) + area(i).theta_s() * (cv.rho_s - cv.rho_g)));
		if (lambda_max < a) lambda_max = a;
	}
}