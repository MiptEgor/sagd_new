#include "energy.h"

energy::energy(mesh &Mesh, const_values CV)
	: area(Mesh)
	, cv(CV)
	, h(cv.h)
	{}

void energy::calc_lay(double t, double tau)
{
//----------------Подсчет граничных значений потоков тепла-------------------------
	if (area.get_left(0).W_l >= 0)
		area.get_left(0).Q_l = cv.c_l * cv.T_bot * area.get_left(0).W_l;
	else
		area.get_left(0).Q_l = cv.c_l * area(0).T * area.get_left(0).W_l;

	if (area.get_left(0).W_g >= 0)
		area.get_left(0).Q_g = cv.c_g * cv.T_bot * area.get_left(0).W_g;
	else
		area.get_left(0).Q_g = cv.c_g * area(0).T * area.get_left(0).W_g;
	
	if (area.get_right(area.get_n() - 1).W_l >= 0)
		area.get_right(area.get_n() - 1).Q_l = cv.c_l * area(area.get_n() - 1).T * area.get_right(area.get_n() - 1).W_l;
	else
		area.get_right(area.get_n() - 1).Q_l = cv.c_l * cv.T_top * area.get_right(area.get_n() - 1).W_l;

	if (area.get_right(area.get_n() - 1).W_g >= 0)
		area.get_right(area.get_n() - 1).Q_l = cv.c_g * area(area.get_n() - 1).T * area.get_right(area.get_n() - 1).W_g;
	else
		area.get_right(area.get_n() - 1).Q_g = cv.c_g * cv.T_top * area.get_right(area.get_n() - 1).W_g;

	area.get_right(area.get_n() - 1).Q_lambda = 0;
	area.get_left(0).Q_lambda = - cv.la * area(0).theta_s() * (area(0).T - cv.T_bot) / h;
	
//--------------Подсчет значений потоков тепла внутри области---------------------------
	for (int i = 0; i < area.get_n() - 1; ++i)
	{
		if (area.get_right(i).W_l >= 0)
			area.get_right(i).Q_l = cv.c_l * area(i).T * area.get_right(i).W_l;
		else
			area.get_right(i).Q_l = cv.c_l * area(i+1).T * area.get_right(i).W_l;

		if (area.get_right(i).W_g >= 0)
			area.get_right(i).Q_g = cv.c_g * area(i).T * area.get_right(i).W_g;
		else
			area.get_right(i).Q_g = cv.c_g * area(i+1).T * area.get_right(i).W_g;

		area.get_right(i).Q_lambda = - cv.la * area(i+1).theta_s() * (area(i+1).T - area(i).T) / h;
	}

//--------------------Расчет значений температуры на новом временном слое
	for (int i = 0; i < area.get_n(); ++i)
	{
		area(i).T = 1. / (cv.c_l * area(i).psi_l + cv.c_s + cv.c_g * area(i).psi_g)  * 
			(area(i).T * (cv.c_l * area(i).old_l + cv.c_s + cv.c_g * area(i).old_g) - tau / h * 
			(area.get_right(i).Q_l + area.get_right(i).Q_g - area.get_left(i).Q_l - area.get_left(i).Q_g
			+ area.get_right(i).Q_lambda - area.get_left(i).Q_lambda));

	}
}