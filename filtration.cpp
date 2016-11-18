	#include "filtration.h"


	filtration::filtration(mesh &Mesh, const_values cv)
		:area(Mesh)
		, cv(cv)
		, h(cv.h)
		, tau(cv.tau)
		{
			find_max();
			std::cout<<max<<" "<<min<<std::endl;
		}

	double filtration::k_perm(double theta)
	{
		return theta * theta;
	}

	double filtration::kappa(double s)
	{
		return k_perm(s) * k_perm((1 - s)) / 
				(k_perm(s) * cv.eta_g + k_perm((1 - s)) * cv.eta_l);

	}
	double filtration::D_kappa(double s)	//Производная kappa по s
	{
		double kl = k_perm(s);
		double kg = k_perm((1 - s));
		double dkl = 2 * s;
		double dkg = 2 * (s - 1);
		return ( (dkl * kg + dkg * kl) / (kl * cv.eta_g + kg * cv.eta_l) 
				- kl * kg * (dkl *cv.eta_g + dkg * cv.eta_l) / pow((kl * cv.eta_g + kg * cv.eta_l), 2) );
	}

	double filtration::lambda(double s)		//Собственное значение (скорость звука)
	{
			return D_kappa(s) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity;
	}

	void filtration::find_max()
	{
		int N = 1e4;
		double hs = 1./N;
		double s = 0;
		max = -4e5;
		min = 4e5;
		for (int i = 0; i < N; ++i)
		{
			if (lambda(s) > max)
				max = lambda(s);
			if (lambda(s) < min)
				min = lambda(s);
			s +=hs;
		}
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
			if (n % 100 == 0)
			{
				area.print_lay(n);
				std::cout<<t<<std::endl;
			}
			calc_lay(t);
			n++;
			t+=tau;
			tau = h / 2 / lambda_max / 2  ;
			if (tau * cv.W_bound_left / h > 0.5) tau = h / cv.W_bound_left / 2 / 2;
		}
	}


	void filtration::calc_lay(double t)
	{
		area.get_left(0).W_l = cv.W_bound_left;
		area.get_right(area.get_n() - 1).W_l = cv.W_bound_right;
		double a = 0;
		lambda_max = 0;
		double min_lambda = 1e10;
		double max_lambda = 0;
		double ds = 1. / area.get_n();
		for (int i = 0; i < area.get_n(); ++i)
		{			//Определение максимального собственного значения (Для условия Куранта)
			a = lambda(i);
			if (max_lambda < a) max_lambda = a;
			if (min_lambda > a) min_lambda = a;
			if (lambda(ds * i) > lambda_max) lambda_max = lambda(ds * i);
		}

		for (int i = 0; i < area.get_n() - 1; ++i)
		{
			double a_r = std::max(lambda(area(i+1).s), lambda(area(i).s));
			double a_l = std::min(lambda(area(i+1).s), lambda(area(i).s));
			area(i).a_r = a_r;
			area(i).a_l = a_l;

			area(i).kappa = kappa(area(i).s) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity;
			if (a_r * a_l >= 0)
			{
				if (a_r > 0) 
					area.get_right(i).W_l = kappa(area(i).s) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity;
				else 
					area.get_right(i).W_l = kappa(area(i+1).s) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity;
			}
			else
			{//-----Вторая часть не обязательна. Почти ничего не меняет
				if(std::abs(area(i).s - area(i+1).s) > 0.4)
					area.get_right(i).W_l = ( (kappa(area(i).s) * max - kappa(area(i+1).s) * min) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity 
							+ max * min * (area(i+1).s - area(i).s) ) / (max - min);
				else
					area.get_right(i).W_l = ( (kappa(area(i).s) * a_r - kappa(area(i+1).s) * a_l) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity 
							+ a_l * a_r * (area(i+1).s - area(i).s) ) / (a_r - a_l);
			
			}
			if (((area(i+1).s == 1.)&(area(i).s == 0.))|((area(i).s == 1.)&(area(i+1).s == 0.)))
				area.get_right(i).W_l = ( (kappa(area(i).s) * max - kappa(area(i+1).s) * min) * cv.K_abs * (cv.rho_l - cv.rho_g) * cv.gravity 
						+ max * min * (area(i+1).s - area(i).s) ) / (max - min);
		}

		for (int i = 0; i < area.get_n(); ++i)
		{
			area(i).s = area(i).s - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau / cv.m;

			//Определение максимального собственного значения
			//a = std::abs(lambda(area(i).s));
			//if (lambda_max < a) lambda_max = a;
		}
	}