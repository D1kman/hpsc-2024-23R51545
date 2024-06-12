#include <cstdio>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;

__global__ void cal(double *u, double *v, double *p, double *b, double *un, double *vn, double *pn, int nx, int ny, int n, int nit, double dx, double dy, double dt, int rho, double nu) {
        int j = blockIdx.x * blockDim.x + threadIdx.x;
	
	printf("%d ", j);
       // for (int n = 0; n < nt; n++) {
	
		if (j % ny != 0 && j % ny != 40) {
        		b[j] = rho*((1/dt) * ((u[j+1] - u[j-1]) / (2*dx) + (v[j+nx] - v[j-nx]) / (2*dy)) - pow(((u[j+1]-u[j-1]) / (2*dx)), 2) -\
                	2*(((u[j+nx] - u[j-nx]) / (2*dy)) * (v[j+1] - v[j-1]) / (2*dx)) - pow(((v[j+nx] - v[j-nx]) / (2*dy)),2));
        	}
		__syncthreads();
	        for (int it = 0; it < nit; it++) {
                	pn[j] = p[j];
                        
                        if (j % ny != 0 && j % ny != 40) {     
                                
                                p[j] = (pow(dy,2) * (pn[j+1] + pn[j-1]) + pow(dx,2) * (pn[j+nx] + pn[j-nx]) - b[j] * pow(dx,2) * pow(dy,2)) / (2*(pow(dx,2) + pow(dy,2)));
                        }
                        
			__syncthreads();
                        if (j < nx) {
                                p[(j*nx)+(ny-1)] = p[(j*nx)+(ny-2)];
                                p[j] = p[ny+j];
                                p[j * ny] = p[(j*ny)+1];
                                p[(ny * (nx-1))+j] = 0;
                        }

                }
		un[j] = u[j];
                vn[j] = v[j];
		__syncthreads();
                        if (j % ny != 0 && j % ny != 40) {

                         
                        
                                u[j] = un[j] - un[j] * dt / dx * (un[j] - un[j-1]) - vn[j] * dt / dy * (un[j] - un[j-ny]) - dt / (2*rho*dx) *\
                                (p[j+1] - p[j-1]) + nu * dt / pow(dx,2) * (un[j+1] - 2 * un[j] + un[j-1]) + nu * dt / pow(dy,2) *\
                                (un[j+ny] - 2 * un[j] + un[j-ny]);
                                v[j] = vn[j] - un[j] * dt / dx * (vn[j] - vn[j-1]) - vn[j] * dt / dy * (vn[j] - vn[j-ny]) - dt / (2*rho*dx) *\
                                (p[j+ny] - p[j-ny]) + nu * dt / pow(dx,2) * (vn[j+1] - 2 * vn[j] + vn[j-1]) + nu * dt / pow(dy,2) *\
                                (vn[j+ny] - 2 * vn[j] + vn[j-ny]);

                        }
		__syncthreads();
                if (j < nx) {
                        u[j] = u[j * ny] = u[j * ny + (ny-1)] = 0.0;
                        u[ny*(nx-1) + j] = 1.0;
                        v[j] = v[ny*(nx-1)+j] = v[j*ny] = v[j*ny + (ny-1)] = 0.0;
                }
 		__syncthreads();


//	}
                
        

}


int main() {
	const int nx = 41;
       	const int ny = 41;
	const int nt = 50;
	const int nit = 50;
	const double dx = 2.0 / (nx - 1);
	const double dy = 2.0 / (ny-1);
	const double dt = 0.01;
	const int rho = 1;
	const double nu = 0.02;
	double* u; cudaMallocManaged(&u, nx*ny*sizeof(double));
        double* v; cudaMallocManaged(&v, nx*ny*sizeof(double));
        double* p; cudaMallocManaged(&p, nx*ny*sizeof(double));
        double* b; cudaMallocManaged(&b, nx*ny*sizeof(double));
        double* un; cudaMallocManaged(&un, nx*ny*sizeof(double));
        double* vn; cudaMallocManaged(&vn, nx*ny*sizeof(double));
        double* pn; cudaMallocManaged(&pn, nx*ny*sizeof(double));
        const int M = 1024;
        const int N = (nx-1)*(ny-1);
	u[nx*ny] = {0};
	v[nx*ny] = {0};
	p[nx*ny] = {0};
	b[nx*ny] = {0};
	pn[nx*ny] = {0};
	un[nx*ny] = {0};
	vn[nx*ny] = {0};

	  ofstream ufile("u.dat");
	  ofstream vfile("v.dat");
	  ofstream pfile("p.dat");


	for (int n = 0; n < nt; n++) {
		printf("walla");
		cal<<<(N+M-1)/M,M>>>(u, v, p, b, un, vn, pn, nx, ny, n, nit, dx, dy, dt, rho, nu);
/*
		for (int j = ny; j < (nx-1)*(ny); j++){
			if (j % ny == 0 || j % ny == 40) {
				
			}
			else {
				b[j] = rho*((1/dt) * ((u[j+1] - u[j-1]) / (2*dx) + (v[j+nx] - v[j-nx]) / (2*dy)) - pow(((u[j+1]-u[j-1]) / (2*dx)), 2) -\
				2*(((u[j+nx] - u[j-nx]) / (2*dy)) * (v[j+1] - v[j-1]) / (2*dx)) - pow(((v[j+nx] - v[j-nx]) / (2*dy)),2));
			}
		}
*/
/*
		for (int it = 0; it < nit; it++) {
			for(int j = 0; j < (nx)*(ny); j++) {
				pn[j] = p[j];
			}

			for(int j = ny; j < (nx-1)*(ny); j++){
				if (j % ny == 0 || j % ny == 40) {
					
				}
				else {
				p[j] = (pow(dy,2) * (pn[j+1] + pn[j-1]) + pow(dx,2) * (pn[j+nx] + pn[j-nx]) - b[j] * pow(dx,2) * pow(dy,2)) / (2*(pow(dx,2) + pow(dy,2)));
				}	
			}

			for (int row = 0; row < nx; row++) {
				p[(row*nx)+(ny-1)] = p[(row*nx)+(ny-2)];
				p[row] = p[ny+row];
				p[row * ny] = p[(row*ny)+1];
				p[(ny * (nx-1))+row] = 0;
			}
		
		}
*/
/*		
		for (int row = 0; row < (nx)*(ny); row++) {
			un[row] = u[row];
			vn[row] = v[row];	
		}
*/
/*		
		for (int j = ny; j < (nx-1)*(ny); j++) {
			if (j % ny == 0 || j % ny == 40) {

                         }
                        else {
				u[j] = un[j] - un[j] * dt / dx * (un[j] - un[j-1]) - vn[j] * dt / dy * (un[j] - un[j-ny]) - dt / (2*rho*dx) *\
				(p[j+1] - p[j-1]) + nu * dt / pow(dx,2) * (un[j+1] - 2 * un[j] + un[j-1]) + nu * dt / pow(dy,2) *\
				(un[j+ny] - 2 * un[j] + un[j-ny]);
				v[j] = vn[j] - un[j] * dt / dx * (vn[j] - vn[j-1]) - vn[j] * dt / dy * (vn[j] - vn[j-ny]) - dt / (2*rho*dx) *\
                                (p[j+ny] - p[j-ny]) + nu * dt / pow(dx,2) * (vn[j+1] - 2 * vn[j] + vn[j-1]) + nu * dt / pow(dy,2) *\
                                (vn[j+ny] - 2 * vn[j] + vn[j-ny]);

			}
		}
*/
/*		
		for(int row = 0; row < nx; row++) {
			u[row] = u[row * ny] = u[row * ny + (ny-1)] = 0.0;
			u[ny*(nx-1) + row] = 1.0;
			v[row] = v[ny*(nx-1)+row] = v[row*ny] = v[row*ny + (ny-1)] = 0.0;
		}
*/	
		
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
 	cudaDeviceSynchronize();
	cudaFree(u);
	cudaFree(v);
	cudaFree(p);
	cudaFree(b);
	cudaFree(un);
	cudaFree(vn);
	cudaFree(pn);  
	}
  ufile.close();
  vfile.close();
  pfile.close();



}

