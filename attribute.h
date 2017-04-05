#ifndef CELLATTR
#define CELLATTR

struct faceAttr
{
	double W_l;
	double W_g;
	double Q_l;
	double Q_g;
	double Q_lambda;
};

struct cellAttr
{
	double T;
	double psi_l;
	double psi_g;
	double eta_l;
	double eta_g;
	double s;
	double a_r;
	double a_l;
	double kappa;
	//double K_perm;
	double x;
	double z;
	double theta_l();
	double theta_s();
	double theta_g();
	double old_l;
	double old_g;
	double lambda_first;
	double lambda_second;
};

#endif //CELLATTR