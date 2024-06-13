#include <cstdio>
#include <math.h>
#include <iostream>
#include <fstream>
#include <cooperative_groups.h>
using namespace cooperative_groups;
using namespace std;

__global__ void calb(double *u, double *v, double *p, double *b, double *un, double *vn, double *pn, int nx, int ny, int nit, double dx, double dy, double dt, double rho, double nu) {
	int j = blockIdx.x;
	int i = threadIdx.x;
	int jp = (j+1)*ny+i;
	int jm = (j-1)*ny+i;
	int ip = j*ny+(i+1);
	int im = j*ny +(i-1);
	int ji = j*ny +i;
	if (j != 0 && j != ny-1 && i != 0 && i != nx-1) {
		b[ji] = rho*((1/dt) * ((u[ip] - u[im]) / (2*dx) + (v[jp] - v[jm]) / (2*dy)) - pow(((u[ip]-u[im]) /\
			(2*dx)), 2) - 2*(((u[jp] - u[jm]) / (2*dy)) * (v[ip] - v[im]) / (2*dx)) - pow(((v[jp] - v[jm]) / (2*dy)),2));
		
	}
}	

__global__ void calp(double *u, double *v, double *p, double *b, double *un, double *vn, double *pn, int nx, int ny, int nit, double dx, double dy, double dt, double rho, double nu) {
        int j = blockIdx.x;
        int i = threadIdx.x;
        int jp = (j+1)*ny+i;
        int jm = (j-1)*ny+i;
        int ip = j*ny+(i+1);
        int im = j*ny +(i-1);
        int ji = j*ny +i;

		pn[ji] = p[ji];
		__syncthreads();
		if (j != 0 && j != ny-1 && i != 0 && i != nx-1) {     

			p[ji] = (pow(dy,2) * (pn[ip] + pn[im]) + pow(dx,2) * (pn[jp] + pn[jm]) - b[ji] * pow(dx,2) * pow(dy,2)) / (2*(pow(dx,2) + pow(dy,2)));
		}
		__syncthreads();
		p[(j*ny)+(ny-1)] = p[(j*ny)+(ny-2)];
		p[i] = p[ny+i];
		p[j*ny] = p[(j*ny)+1];
		p[(ny * (nx-1))+i] = 0.0;

		
}	
__global__ void caluv(double *u, double *v, double *p, double *b, double *un, double *vn, double *pn, int nx, int ny, int nit, double dx, double dy, double dt, double rho, double nu) {
        int j = blockIdx.x;
        int i = threadIdx.x;
        int jp = (j+1)*ny+i;
        int jm = (j-1)*ny+i;
        int ip = j*ny+(i+1);
        int im = j*ny +(i-1);
        int ji = j*ny +i;
	un[ji] = u[ji];
	vn[ji] = v[ji];
	if (j != 0 && j != ny-1 && i != 0 && i != nx-1) {


		u[ji] = un[ji] - un[ji] * dt / dx * (un[ji] - un[im]) - vn[ji] * dt / dy * (un[ji] - un[jm]) - dt / (2*rho*dx) *\
			(p[ip] - p[im]) + nu * dt / pow(dx,2) * (un[ip] - 2 * un[ji] + un[im]) + nu * dt / pow(dy,2) *\
			(un[jp] - 2 * un[ji] + un[jm]);
		v[ji] = vn[ji] - un[ji] * dt / dx * (vn[ji] - vn[im]) - vn[ji] * dt / dy * (vn[ji] - vn[jm]) - dt / (2*rho*dx) *\
			(p[jp] - p[jm]) + nu * dt / pow(dx,2) * (vn[ip] - 2 * vn[ji] + vn[im]) + nu * dt / pow(dy,2) *\
			(vn[jp] - 2 * vn[jp] + vn[jm]);

	}
	__syncthreads();
	u[i] = u[j*ny] = u[j*ny +(ny-1)] = 0.0;
	u[ny*(nx-1) + i] = 1.0;
	v[i] = v[ny*(nx-1)+i] = v[j*ny] = v[j*ny + (ny-1)] = 0.0;

}


int main() {
	const int nx = 41;
	const int ny = 41;
	const int nt = 120;
	const int nit = 50;
	const double dx = 2.0 / (nx - 1);
	const double dy = 2.0 / (ny-1);
	const double dt = 0.01;
	const double rho = 1.0;
	const double nu = 0.02;
	double* u; cudaMallocManaged(&u, nx*ny*sizeof(double));
	double* v; cudaMallocManaged(&v, nx*ny*sizeof(double));
	double* p; cudaMallocManaged(&p, nx*ny*sizeof(double));
	double* b; cudaMallocManaged(&b, nx*ny*sizeof(double));
	double* un; cudaMallocManaged(&un, nx*ny*sizeof(double));
	double* vn; cudaMallocManaged(&vn, nx*ny*sizeof(double));
	double* pn; cudaMallocManaged(&pn, nx*ny*sizeof(double));

	cudaMemset(u, 0, nx*ny*sizeof(double));
	cudaMemset(v, 0, nx*ny*sizeof(double));
	cudaMemset(p, 0, nx*ny*sizeof(double));
	cudaMemset(b, 0, nx*ny*sizeof(double));
	cudaMemset(un, 0, nx*ny*sizeof(double));
	cudaMemset(vn, 0, nx*ny*sizeof(double));
	cudaMemset(pn, 0, nx*ny*sizeof(double));

	ofstream ufile("u.dat");
	ofstream vfile("v.dat");
	ofstream pfile("p.dat");

	for (int n = 0; n < nt; n++) {
		calb<<<nx,ny>>>(u, v, p, b, un, vn, pn, nx, ny, nit, dx, dy, dt, rho, nu);

		cudaDeviceSynchronize();

		for (int it = 0; it < nit; it++) {
			calp<<<nx, ny>>>(u, v, p, b, un, vn, pn, nx, ny, nit, dx, dy, dt, rho, nu);

			cudaDeviceSynchronize();
		}

		caluv<<<nx, ny>>>(u, v, p, b, un, vn, pn, nx, ny, nit, dx, dy, dt, rho, nu);

		cudaDeviceSynchronize();
		if (n % 10 == 0) {
			for (int j=0; j<ny; j++)
				for (int i=0; i<nx; i++)
					ufile << u[j*ny+i] << " ";
			ufile << "\n";
			for (int j=0; j<ny; j++)
				for (int i=0; i<nx; i++)
					vfile << v[(j*ny)+i] << " ";
			vfile << "\n";
			for (int j=0; j<ny; j++)
				for (int i=0; i<nx; i++)
					pfile << p[j*ny+i] << " ";
			pfile << "\n";
		}

	}
	cudaFree(u);
	cudaFree(v);
	cudaFree(p);
	cudaFree(b);
	cudaFree(un);
	cudaFree(vn);
	cudaFree(pn); 
	ufile.close();
	vfile.close();
	pfile.close();

}

