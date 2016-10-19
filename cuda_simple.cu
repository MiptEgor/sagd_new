#include <cmath>
#include <iostream>
#include <fstream>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <algorithm>
#include <thrust/extrema.h>
#define BLOCK_SIZE 16

float f(float g, float l)
{
	return (962361*(pow(-74500000*(-2 + g)*g + 220*l + 
     6*(31333357 - 3250013*g)*g*l + (-110 + 39000039*g)*l*l,2) / 2500 -
   1600*g*l*(1 + g + l)*(5811*g*g - 110*(-2 + l)*(149 + 39*l) - 
      2*g*(14006 + 1521*l))))/(100 * pow((1 + g + l),8));
}
__device__ float k(float tetta)
{
	return tetta * tetta;
}

__device__ float k_shtrih(float tetta)
{
	return tetta * tetta;
}

__device__ float W1(float l, float g, float rg, float rl, float xi) 
{
	return   l * (l - 2 - 2 * g) * (g * (rg - rl) - rl + 1) / pow(1 + g + l, 4);
}
__device__ float W2(float l, float g, float rg, float rl, float xi)
{
	return  - l * l * ((1 - 2 * g + l) * rg + (2 + 2  * g - l) * rl - 3) / pow(1 + g + l, 4);
}
__device__ float W3(float l, float g, float rg, float rl, float xi)
{
	return  g * g * ((g - 2 - 2 * l) * rg - (1 + g - 2 * l) * rl + 3) / pow(1 + g + l, 4) / xi;
}
__device__ float W4(float l, float g, float rg, float rl, float xi)
{
	return + g * (g - 2 * l - 2) * (rg + rg * l - l * rl - 1) / pow(1 + g + l, 4) / xi;
} 

__global__ void parallel_calculate(thrust::device_ptr<float> huge_dev, int N, float rg, float rl, float xi, float h)
{
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int j = blockIdx.y * blockDim.y + threadIdx.y;
	float g = i * h;
	float l = j * h;
	//float b = W1(l, g , rg, rl, xi) + W4(l, g , rg, rl, xi);
	if ((i < N)&&(j < N))
	{
		float D = powf( W1(l, g , rg, rl, xi) - W4(l, g , rg, rl, xi), 2) + 4 * W3(l, g , rg, rl, xi) * W2(l, g , rg, rl, xi) 
		 - 4 * W1(l, g , rg, rl, xi) * W4(l, g , rg, rl, xi);
	//float dd = powf(fabs(D), 0.5f);
		huge_dev[i + j * N] = D;
		if (D>0)
			huge_dev[i + j * N] = 0;
	/*if (D < dev[2])
	{
		dev[0] = g;
		dev[1] = l;
		dev[2] = D;
		dev[3] = b;
		dev[4] = dd;
	}*/
	}
	__syncthreads();

}


template<typename T>
void put_be(std::ofstream &f, const T val) {
    union {
        T value;
        char bytes[sizeof(T)];
    } x;
    x.value = val;
    std::reverse(x.bytes, x.bytes + sizeof(T));
    f.write(x.bytes, sizeof(T));
}

void print(int nx, int ny, thrust::host_vector<float>::iterator huge_host, float h)
	{

		int Nx = nx;
		int Ny = ny;
	    char path[1024];
	    sprintf(path, "res.vtk");
	    std::ofstream f(path, std::ios::binary);

	    f << "# vtk DataFile Version 3.0" << std::endl;
	    f << "Comment" << std::endl;
	    f << "BINARY" << std::endl;
	    f << "DATASET RECTILINEAR_GRID" << std::endl;
	    f << "DIMENSIONS " << Nx << " " << Ny << " 1" << std::endl;
	    f << "X_COORDINATES " << Nx << " float" << std::endl;
	    for (size_t i = 0; i < Nx; i++)
	        put_be<float>(f, h * i);
	    f << "Y_COORDINATES " << Ny << " float" << std::endl;
	    for (size_t j = 0; j < Ny; j++)
	        put_be<float>(f, h * j);
	    f << "Z_COORDINATES 1 float" << std::endl;
	    put_be<float>(f, 0);

	    f << "POINT_DATA " << Nx * Ny << std::endl;
	    f << "SCALARS p float\nLOOKUP_TABLE default" << std::endl;
	    for (size_t j = 0; j < Ny; j++)
	        for (size_t i = 0; i < Nx; i++)
	            put_be<float>(f, huge_host[i + j * Nx]);
	    f.close();

	}

int main(int argc, char const *argv[])
{
	std::cout.precision(10);
	int N = 2e3;
	float h = 1.e2 / N;
	dim3 block(BLOCK_SIZE, BLOCK_SIZE);
	dim3 nblocks((int)(N) / BLOCK_SIZE + 1, (int)(N) / BLOCK_SIZE + 1);
	thrust::device_vector<float> dev(5, 0.f);
	thrust::host_vector<float> host(5, 0.f);
	thrust::device_vector<float> huge_dev(N * N, 0.f);
	thrust::host_vector<float> huge_host(N * N, 0.f);
	float rl = 0.3;
	float rg = 0.02;
	float xi = 1.e-4;
	std::cout<<"g"<<" "<<"l"<<" "<<"D"<<" "<<"b"<<" "<<"dd" <<std::endl;
	parallel_calculate<<<nblocks, block>>>(huge_dev.data(), N, rg, rl, xi, h);
	thrust::copy(huge_dev.begin(), huge_dev.end(), huge_host.begin());
	//std::cout<<huge_host[0]<<" "<<huge_host[1]<<" "<<huge_host[2]<<" "<<huge_host[3]<<" "<<huge_host[4] <<std::endl;
	print(N, N, huge_host.data(), h);
	float min = *(thrust::min_element(huge_dev.begin(), huge_dev.end()));
	std::cout<<min<<std::endl;
	/*xi = 1e-2;
	parallel_calculate<<<nblocks, block>>>(dev.data(), N, rg, rl, xi, h);
	thrust::copy(dev.begin(), dev.end(), host.begin());
	std::cout<<host[0]<<" "<<host[1]<<" "<<host[2]<<" "<<host[3]<<" "<<host[4] <<std::endl;
	*/
	std::cout<<"end"<<std::endl;
	return 0;
}