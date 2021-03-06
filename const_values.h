#ifndef CONST_VALUES
#define CONST_VALUES

struct const_values
{


const_values()
	: T(2.) 
	, T0(300.0)					//Температура среды и верхней стенки
	, T1(500.0)					//Температура нижней стенки		(К)
	, length(10.0)  			//Длина области в метрах		(м)
	, time_const ( 1e10)			//Время эксперимента 			(с)
	, tau(5e0)	//Величина шага по времени		(с) 
	, h(0.01)				//Величина шага по расстоянию	(м)
	, eta_l(1e-3)
	//, eta_g(18.3e-5 )					//Коэффициент Вязости			(Па*с)
	//, eta_g(5e-4)
	, eta_g(1e-3)
	, la(0.6)					//Коэффициент Теплопроводности	(Джpar/(м*с*К))
	, rho_l(1000.0) 				//Плотность жидкости			(Кг/м3)
	, rho_s(1400.0) 				//Плотность скелета				(Кг/м3)
	, rho_g(500)
	, c_l(1e3)				//Теплоемкость флюида			(Дж/К)
	, c_s(2e3)				//Теплоемкость скелета
	, gravity(-9.8)				//Постоянная свободного падения	(м/с2)
	, K_abs(1e-11)				//Абсолютная проницаемость
	, psi0(2.0)					//Насыщенность на границе
	, psi_crit(0.1/0.9)			//Критическое значение насыщенности (относительной)
	, alpha(0.02)				//Коэффициент в показателе экспоненты для расчета вязкости
	, Tcrit(400)				//Температура в экспоненте для расчечта вязкоти
	, m(1.)					//Пористость
	, s(0.5)					//Начальная насыщенность жидкостью
	, W_bound_left(-5e-7)
	, W_bound_right(0)
	{

	}
	double T;
	double T0;					//Температура среды и верхней стенки
	double T1;					//Температура нижней стенки		(К)
	double length;  		//Длина области в метрах		(м)
	double time_const;		//Время эксперимента 			(с)
	double tau;	//Величина шага по времени		(с) 
	double h;				//Величина шага по расстоянию	(м)
	double eta_l;
	double eta_g;
	//double eta; 			//Коэффициент Вязости			(Па*с)
	double la;				//Коэффициент Теплопроводности	(Дж/(м*с*К))
	double rho_l;			//Плотность жидкости			(Кг/м3)
	double rho_s; 			//Плотность скелета				(Кг/м3)
	double rho_g;
	double c_l;				//Теплоемкость	флюида				(Дж/К)
	double c_s;				//Теплоемкость скелета				(Дж/К)
	double gravity;			//Постоянная свободного падения	(м/с2)
	double K_abs;			//Абсолютная проницаемость
	double psi0;			//Насыщенность на границе
	double psi_crit;		//Критическое значение насыщенности (относительной)
	double alpha;		//Коэффициент в показателе экспоненты для расчета вязкости
	double Tcrit;		//Температура в экспоненте для расчечта вязкоти
	double m;
	double s;
	double W_bound_left;
	double W_bound_right;

};


#endif