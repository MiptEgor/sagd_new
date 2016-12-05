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
		cell[0].eta_l = cv.eta_l;
		cell[0].eta_g = cv.eta_g;
		cell[0].psi_l = cv.psi_l;
		cell[0].psi_g = cv.psi_g;
		for (int i = 1; i < nx; ++i)
		{
			cell[i].x = cell[i - 1].x + dx;
			cell[i].eta_l = cv.eta_l;
			cell[i].eta_g = cv.eta_g;
			cell[i].psi_l = cv.psi_l;
			cell[i].psi_g = cv.psi_g;
		}
	}

void mesh::print_lay(int n)
{
	std::fstream fs;
	char buf[128];
	sprintf(buf,"./out/output%07d.csv",n);
	fs.open(buf, std::fstream::out);
	fs << "x, z, theta_l, theta_g, a_r, a_l, W_l, W_g" << std::endl;
	cell[0].z = 0;
	for (int i = 1; i < nx; ++i)
	{
		//Вычисление Эйлеровой коородинаты
		//Почему не вместе с выводом - на случай если не каждую ячейку выводить
		cell[i].z = cell[i-1].z + dx * (1./cell[i].theta_s());
	}

	for (int i = 0; i < nx; i+=2)
	{
		fs<<cell[i].x<<"," << cell[i].z << "," <<cell[i].theta_l() <<"," <<cell[i].theta_g()<<"," 
		<< cell[i].a_r<<"," << cell[i].a_l<<"," <<face[i].W_l<<"," <<face[i].W_g<<"\n";
	}
	fs.close();	
}