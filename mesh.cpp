#include "attribute.cpp"
#include "mesh.h"
#include "const_values.h"
mesh::mesh(const_values cv)
	: X(cv.length)
	, dx(cv.h)
	, nx(X/dx + 1)
	, cell(nx)
	, face(nx + 1)
	{
		cell[0].x = 0;
		for (int i = 1; i < nx; ++i)
		{
			cell[i].x = cell[i - 1].x + dx;
			cell[i].eta_l = 1e-3;
			cell[i].eta_g = 18.3e-6;
		}

	}
void mesh::print_lay(int n)
{
	std::fstream fs;
	char buf[128];
	sprintf(buf,"./out/output%07d.csv",n);
	fs.open(buf, std::fstream::out);
	fs << "x, eta_l, eta_g" << std::endl;
	for (int i = 0; i < nx; ++i)
	{
		fs<<cell[i].x <<"," <<cell[i].theta_l() <<"," <<cell[i].theta_g() <<std::endl;
	}
	fs.close();	
}