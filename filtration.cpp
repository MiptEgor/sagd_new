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

	double filtration::kappa(int i)
	{
		return cv.K_abs * k_perm(area(i).s) * k_perm((1 - area(i).s)) / 
				(k_perm(area(i).s) * area(i).eta_g + k_perm((1 - area(i).s)) * area(i).eta_l);

	}

	double filtration::lambda(int i)
	{
			double kl = k_perm(area(i).s);
			double kg = k_perm((1 - area(i).s));
			double dkl = 2 * area(i).s;
			double dkg = 2 * (area(i).s - 1);
			return ( kappa(i) * (dkl / kl + dkg / kg) - kappa(i) / (kl * cv.eta_g + kg * cv.eta_l) * (dkl * cv.eta_g + dkg * cv.eta_l) ) * (cv.rho_l - cv.rho_g) * cv.gravity;
			//cv.m * cv.m *
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
			tau = h / 2 / lambda_max / 2 / 10;
			if (tau * cv.W_bound_left / h > 0.5) tau = h / cv.W_bound_left / 2;
		}
	}


	void filtration::calc_lay(double t)
	{
		area.get_left(0).W_l = cv.W_bound_left;
		area.get_right(area.get_n() - 1).W_l = cv.W_bound_right;
		double a = 0;
		lambda_max = 0;
		double s_ = pow(cv.eta_l/ cv.eta_g, 1./3) / (1 + pow(cv.eta_l/ cv.eta_g, 1./3));
		for (int i = 0; i < area.get_n() - 1; ++i)
		{
			
			if ((area(i).s <= s_)&(area(i+1).s <= s_))
				area.get_right(i).W_l = kappa(i) * (cv.rho_l - cv.rho_g) * cv.gravity;
			if ((area(i).s >= s_)&(area(i+1).s >= s_))
				area.get_right(i).W_l = kappa(i+1) * (cv.rho_l - cv.rho_g) * cv.gravity;
			
			if ((area(i).s < s_)&(area(i+1).s > s_))
				area.get_right(i).W_l = cv.K_abs * k_perm(s_ * cv.m) * k_perm((1 - s_) * cv.m) / 
				(k_perm(s_ * cv.m) * area(i).eta_g + k_perm((1 - s_) * cv.m) * area(i).eta_l) * (cv.rho_l - cv.rho_g) * cv.gravity;

			if ((area(i).s < s_)&(area(i+1).s > s_))
			{
				if (kappa(i) < kappa(i + 1))
					area.get_right(i).W_l = kappa(i) * (cv.rho_l - cv.rho_g) * cv.gravity;
				else
					area.get_right(i).W_l = kappa(i + 1) * (cv.rho_l - cv.rho_g) * cv.gravity;
			}

		}

		for (int i = 0; i < area.get_n(); ++i)
		{
			area(i).s = area(i).s - (area.get_right(i).W_l - area.get_left(i).W_l) / h * tau / cv.m;

			//Определение максимального собственного значения
			a = std::abs(lambda(i));
			if (lambda_max < a) lambda_max = a;
		}
	}