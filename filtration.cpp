#include "filtration.h"


filtration::filtration(mesh &Mesh, const_values cv)
:area(Mesh)
, cv(cv)
, h(cv.h)
{
	find_max();
	std::cout << min << " " << max<<std::endl;
}
double filtration::calc_tau()
{
	return h / 2 / lambda_max / 4;
}

double filtration::k_perm(double theta)
{
	return theta * theta;
}

void filtration::initial()
{
	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).psi_l = cv.psi_l;
		area(i).psi_g = cv.psi_g;
	}

}
void filtration::find_max()
{
	int N = 1000;
	max = -1e3;
	min = 1e3;
	for (int i = 0; i < N; ++i)
	{
		for (int j = 0; j < N-i; ++j)
		{
			if(i != j)
			{
				double theta_l = 1. * i / N;
				double theta_g = 1. * j / N;
				double theta_s = 1. - theta_l - theta_g;
				double psi_l = theta_l / theta_s;
				double psi_g = theta_g / theta_s;
				if (first_lambda(psi_l, psi_g) < min)
					min = first_lambda(psi_l, psi_g);
				if (second_lambda(psi_l, psi_g) > max)
					max = second_lambda(psi_l, psi_g);
			}	
		}
	}
}

//------------------Функции для определения частныъ производных W_l и W_G по psi_l и psi_g
double filtration::W_ll(double l, double g)
{

	return - cv.K_abs * l * (2 * g + 2 - l) *(g * (cv.rho_g - cv.rho_l) - cv.rho_l + cv.rho_s) * cv.gravity / cv.eta_l / pow(g + l + 1, 4); 
}

double filtration::W_lg(double l, double g)
{
	return - cv.K_abs * l * l * (cv.rho_g * (l + 1 - 2 * g) + cv.rho_l * (2 * g + 2 - l) - cv.rho_s) * cv.gravity / cv.eta_l / pow(g + l + 1, 4);
}

double filtration::W_gl(double l, double g)
{
	return - cv.K_abs * g * g * (cv.rho_g * (2 * l + 2 - g) + cv.rho_l * (g - 2 * l + 1) - cv.rho_s) * cv.gravity / cv.eta_g / pow(g + l + 1, 4);
}

double filtration::W_gg(double l, double g)
{
	return -cv.K_abs * g * (2 * l + 2 - g) * (l * (cv.rho_l - cv.rho_g) - cv.rho_g + cv.rho_s) * cv.gravity / cv.eta_g / pow(g + l + 1, 4);
}

double filtration::first_lambda(double l, double g)
{	
	double D = pow(W_ll(l, g) + W_gg(l, g), 2) - 4 * (W_ll(l, g) * W_gg(l, g) - W_lg(l, g) * W_gl(l, g));
	return (W_ll(l, g) + W_gg(l, g) - std::sqrt(D)) / 2;	
}

double filtration::second_lambda(double l, double g)
{	
	double D = pow(W_ll(l, g) + W_gg(l, g), 2) - 4 * (W_ll(l, g) * W_gg(l, g) - W_lg(l, g) * W_gl(l, g));
	return (W_ll(l, g) + W_gg(l, g) + std::sqrt(D)) / 2;	
}

double filtration::W_l(double l, double g)
{
	return - l * l / cv.eta_l / pow(1 + l + g, 3) * cv.K_abs * ((cv.rho_s - cv.rho_l) + g * (cv.rho_g - cv.rho_l)) * cv.gravity;
}

double filtration::W_g(double l, double g)
{
	return - g * g / cv.eta_g / pow(1 + l + g, 3) * cv.K_abs * ((cv.rho_s - cv.rho_g) + l * (cv.rho_l - cv.rho_g)) * cv.gravity;
}


void filtration::calc_lay(double t, double tau)
{

	area.get_left(0).W_l = cv.W_bound_left_L;
	area.get_left(0).W_g = cv.W_bound_left_G;
	area.get_right(area.get_n() - 1).W_l = cv.W_bound_right_L;
	area.get_right(area.get_n() - 1).W_g = cv.W_bound_right_G;
	
//----------Определение максимального собственного значения (Для условия Куранта)----------------------
	lambda_max = 0;
	for (int i = 0; i < area.get_n(); ++i)
	{
		if (std::abs(first_lambda(area(i).psi_l, area(i).psi_g)) > lambda_max) lambda_max = std::abs(first_lambda(area(i).psi_l, area(i).psi_g));
		if (std::abs(second_lambda(area(i).psi_l, area(i).psi_g)) > lambda_max) lambda_max = std::abs(second_lambda(area(i).psi_l, area(i).psi_g));
	}

	for (int i = 0; i < area.get_n() - 1; ++i)
	{
		//double a_r = std::max(second_lambda(area(i).psi_l, area(i).psi_g), second_lambda(area(i+1).psi_l, area(i+1).psi_g));
		//double a_l = std::min(first_lambda(area(i).psi_l, area(i).psi_g), first_lambda(area(i+1).psi_l, area(i+1).psi_g));
		double a_r = max;
		double a_l = min;
		area(i).a_r = a_r;
		area(i).a_l = a_l;
		
		if (a_r * a_l > 0)
		{
			if (a_r > 0) 
			{
				area.get_right(i).W_l = W_l(area(i).psi_l, area(i).psi_g);
				area.get_right(i).W_g = W_g(area(i).psi_l, area(i).psi_g);
			}
			else 
			{
				area.get_right(i).W_l = W_l(area(i+1).psi_l, area(i+1).psi_g);
				area.get_right(i).W_g = W_g(area(i+1).psi_l, area(i+1).psi_g);
			}
		
		}
		else
		{
			area.get_right(i).W_l =( W_l(area(i).psi_l, area(i).psi_g) * max - W_l(area(i+1).psi_l, area(i+1).psi_g) * min + 
				max * min * (area(i+1).psi_l - area(i).psi_l) ) / (max - min);

			area.get_right(i).W_g =( W_g(area(i).psi_l, area(i).psi_g) * max - W_g(area(i+1).psi_l, area(i+1).psi_g) * min + 
				max * min * (area(i+1).psi_g - area(i).psi_g) ) / (max - min);
		}
	}

	for (int i = 0; i < area.get_n(); ++i)
	{
		
		area(i).psi_l = area(i).psi_l - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau;
		
		area(i).psi_g = area(i).psi_g - (area.get_right(i).W_g - area.get_left(i).W_g) / h * tau;
	}
}