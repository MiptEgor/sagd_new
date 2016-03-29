#ifndef MESH
#define MESH
#include "const_values.h"
#include <vector>
#include "attribute.h"
#include <stdio.h>
#include <fstream>
class mesh
{
private:
	double X, dx;
	int nx;
	//std::vector<double> x;
	std::vector<cellAttr> cell;
	std::vector<faceAttr> face;
public:
	mesh(const_values cv);
	void print_lay(int n);

	const int & get_n()
	{
		return nx;
	}

	faceAttr & get_left(int i)
	{
		return face[i];
	}
	faceAttr & get_right(int i)
	{
		return face[i + 1];
	}
	cellAttr & operator()(int i)
	{
		return cell[i];
	} 
};

#endif //MESH