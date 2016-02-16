#ifndef CELLATTR
#define CELLATTR

struct faceAttr
{
	double W_l;
	double W_g;
};

struct cellAttr
{
	double T;
	double psi_l;
	double psi_g;
	double eta_l;
	double eta_g;
	//double K_perm;
	double x;
	double theta_l();
	double theta_s();
	double theta_g();
};

#endif //CELLATTR