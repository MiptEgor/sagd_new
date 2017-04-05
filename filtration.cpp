#include "filtration.h"


filtration::filtration(mesh &Mesh, const_values cv)
:area(Mesh)
, cv(cv)
, h(cv.h)
{
	find_max();
	std::cout << min << " " << max<<std::endl;
}

double filtration::viscosity(double T)
{
	return std::exp(-cv.alpha * (T - cv.Tcrit));
}
double filtration::calc_tau()
{
	return h / 2 / lambda_max * 4 ;
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
		{// до N - i чтобы theta_s не уходила ниже нуля
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
//------------------Функции для определения частныъ производных W_l и W_g по psi_l и psi_g
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
	//assert(D>=0);
	return (W_ll(l, g) + W_gg(l, g) - std::sqrt(D)) / 2;	
}

double filtration::second_lambda(double l, double g)
{	
	double D = pow(W_ll(l, g) + W_gg(l, g), 2) - 4 * (W_ll(l, g) * W_gg(l, g) - W_lg(l, g) * W_gl(l, g));
	//assert (D>=0);
	return (W_ll(l, g) + W_gg(l, g) + std::sqrt(D)) / 2;	
}

void filtration::test_D(double l, double g, int i, double t)
{
	double D = pow(W_ll(l, g) + W_gg(l, g), 2) - 4 * (W_ll(l, g) * W_gg(l, g) - W_lg(l, g) * W_gl(l, g));
	if ((D < 0))
	{
		std::cout<<"D<0, t = " <<t<<" i = "<<i<<" "<<std::endl;
		std::terminate();
	}
}

double filtration::W_l(double l, double g, double eta)
{
	return - l * l / eta / pow(1 + l + g, 3) * cv.K_abs * ((cv.rho_s - cv.rho_l) + g * (cv.rho_g - cv.rho_l)) * cv.gravity;
}

double filtration::W_g(double l, double g, double eta)
{
	return - g * g / eta / pow(1 + l + g, 3) * cv.K_abs * ((cv.rho_s - cv.rho_g) + l * (cv.rho_l - cv.rho_g)) * cv.gravity;
}

void filtration::calc_lay(double t, double tau)
{
	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).eta_l = cv.eta_l;
		area(i).eta_g = cv.eta_g;
		//area(i).eta_l = 0.1 * viscosity(area(i).T);
		//area(i).eta_g = 0.01 * viscosity(area(i).T);
	}

	area.get_left(0).W_l = cv.W_bound_left_L;
	area.get_left(0).W_g = cv.W_bound_left_G;
	area.get_right(area.get_n() - 1).W_l = cv.W_bound_right_L;
	area.get_right(area.get_n() - 1).W_g = cv.W_bound_right_G;
	
//----------Определение максимального собственного значения (Для условия Куранта)----------------------
	lambda_max = std::max(std::abs(max), std::abs(min));
	/*
	lambda_max = 0;
	for (int i = 0; i < area.get_n(); ++i)
	{
		if (std::abs(first_lambda(area(i).psi_l, area(i).psi_g)) > lambda_max) lambda_max = std::abs(first_lambda(area(i).psi_l, area(i).psi_g));
		if (std::abs(second_lambda(area(i).psi_l, area(i).psi_g)) > lambda_max) lambda_max = std::abs(second_lambda(area(i).psi_l, area(i).psi_g));
	}
	*/
	

	area.update();
	for (int i = 0; i < area.get_n() - 1; ++i)
	{
		//double a_r = std::max(2 * second_lambda(area(i).old_l, area(i).old_g), 2 * second_lambda(area(i+1).old_l, area(i+1).old_g));
		//double a_l = std::min(2 * first_lambda(area(i).old_l, area(i).old_g), 2 * first_lambda(area(i+1).old_l, area(i+1).old_g));
		//a_l = std::min (a_l, 0.5 * first_lambda(area(i).old_l, area(i).old_g));
		//a_l = std::min (a_l, 0.5 * first_lambda(area(i+1).old_l, area(i+1).old_g));
		//a_r = std::max (a_r, 0.5 * second_lambda(area(i).old_l, area(i).old_g));
		//a_r = std::max (a_r, 0.5 * second_lambda(area(i+1).old_l, area(i+1).old_g));
		double a_r = max;
		double a_l = min;
		area(i).a_r = a_r;
		area(i).a_l = a_l;

		area(i).lambda_first = first_lambda(area(i).psi_l, area(i).psi_g);
		area(i).lambda_second = second_lambda(area(i).psi_l, area(i).psi_g);
		//if (area(i).lambda_second !=area(i).lambda_second) std::cout<<"lambda = nan, t = " <<t<<" i = "<<i<<" "<<std::endl;
		test_D(area(i).psi_l, area(i).psi_g, i, t);
		//---------- Если собственные значения одного знака то сносим справа или слева, в зависимости от знака
		if (a_r * a_l > 0)
		{
			if (a_r > 0) 
			{
				area.get_right(i).W_l = W_l(area(i).old_l, area(i).old_g, area(i).eta_l);
				area.get_right(i).W_g = W_g(area(i).old_l, area(i).old_g, area(i).eta_g);
			}
			else 
			{
				area.get_right(i).W_l = W_l(area(i+1).old_l, area(i+1).old_g, area(i+1).eta_l);
				area.get_right(i).W_g = W_g(area(i+1).old_l, area(i+1).old_g, area(i+1).eta_g);
			}
		}
		else
		{
			area.get_right(i).W_l =( W_l(area(i).old_l, area(i).old_g, area(i).eta_l) * a_r - W_l(area(i+1).old_l, area(i+1).old_g, area(i+1).eta_l) * a_l + 
				a_r * a_l * (area(i+1).old_l - area(i).old_l) ) / (a_r - a_l);

			area.get_right(i).W_g =( W_g(area(i).old_l, area(i).old_g, area(i).eta_g) * a_r - W_g(area(i+1).old_l, area(i+1).old_g, area(i+1).eta_g) * a_l + 
				a_r * a_l * (area(i+1).old_g - area(i).old_g) ) / (a_r - a_l);
		}
	
		area(i).psi_l = area(i).psi_l - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau;
		
		area(i).psi_g = area(i).psi_g - (area.get_right(i).W_g - area.get_left(i).W_g) / h * tau;

		//----------------Проверка на предельность объемной доли твердой фазы--------------------
		
		if (area(i).theta_s() > cv.theta_crit)
		{
			double l = area(i).psi_l;
			double g = area(i).psi_g;
			double l_old = area(i).old_l;
			double g_old = area(i).old_g;
			double c = (1. / cv.theta_crit  - 1 - l_old - g_old) / (l - l_old + g - g_old);
			area(i).psi_l = l_old + c * (l - l_old);
			area(i).psi_g = g_old + c * (g - g_old);
			area.get_right(i).W_l = area.get_left(i).W_l + (area(i).old_l - area(i).psi_l) * h / tau;
			area.get_right(i).W_g = area.get_left(i).W_g + (area(i).old_g - area(i).psi_g) * h / tau;
			//std::cout<<"EPTA"<<std::endl;
		}
		

	}
	
	area(area.get_n() - 1).psi_l = area(area.get_n() - 1).psi_l - (area.get_right(area.get_n() - 1).W_l * 0. - area.get_left(area.get_n() - 1).W_l) / h * tau;
		
	area(area.get_n() - 1).psi_g = area(area.get_n() - 1).psi_g - (area.get_right(area.get_n() - 1).W_g * 0. - area.get_left(area.get_n() - 1).W_g) / h * tau;


	/*
	for (int i = 0; i < area.get_n(); ++i)
	{
		
		area(i).psi_l = area(i).psi_l - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau;
		
		area(i).psi_g = area(i).psi_g - (area.get_right(i).W_g - area.get_left(i).W_g) / h * tau;
	}
	*/
}