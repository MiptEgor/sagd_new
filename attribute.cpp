#include "attribute.h"
double cellAttr::theta_s()
{
	return 1 /(1 + psi_l + psi_g); 
}

double cellAttr::theta_l()
{
	//return theta_s() * psi_l;
	return s;
}

double cellAttr::theta_g()
{
	//return theta_s() * psi_g;
	return 1 - s;
}