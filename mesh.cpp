#include "mesh.h"
#include "const_values.h"
mesh::mesh(const_values cv)
	: X(cv.length)
	, dx(cv.h)
	, nx(X/dx + 1)
	, cell(nx)
	, face(nx + 1)
	{
		double s_top = 1.;
		double s_bot = 0.;
		cell[0].x = 0;
		cell[0].eta_l = cv.eta_l;
		cell[0].eta_g = cv.eta_g;
		cell[0].s = s_bot;
		for (int i = 1; i < nx; ++i)
		{
			cell[i].x = cell[i - 1].x + dx;
			cell[i].eta_l = cv.eta_l;
			cell[i].eta_g = cv.eta_g;
			cell[i].s = s_top;
		}
		for (int i = 1; i < nx/2; ++i)
		{
			cell[i].s = s_bot;
		}
	}

void mesh::print_lay(int n)
{
	std::fstream fs;
	char buf[128];
	sprintf(buf,"./out/output%07d.csv",n);
	fs.open(buf, std::fstream::out);
	fs << "x, s, W_l, a_r, a_l, kappa" << std::endl;
	//double z = 0;
	for (int i = 0; i < nx; i+=5)
	{
		//z = z + dx * (1./cell[i].theta_s()); //Эйлерова координата
		fs<<cell[i].x<<"," <<cell[i].s <<"," << face[i].W_l <<"," << cell[i].a_r<<"," << cell[i].a_l<<"," <<cell[i].kappa<<"\n";
	}
	fs.close();	
}