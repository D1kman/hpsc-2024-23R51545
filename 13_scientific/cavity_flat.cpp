#include <cstdio>
#include <math.h>
#include <iostream>
#include <fstream>
using namespace std;
int main() {
	const int nx = 41;
       	const int ny = 41;
	const int nt = 5;
	const int nit = 50;
	const double dx = 2.0 / (nx - 1);
	const double dy = 2.0 / (ny-1);
	const double dt = 0.01;
	const int rho = 1;
	const double nu = 0.02;
	double x[nx];
	double y[ny];
	double u[nx*ny] = {0};
	double v[nx*ny] = {0};
	double p[nx*ny] = {0};
	double b[nx*ny] = {0};
	double pn[nx*ny] = {0};
	double un[nx*ny] = {0};
	double vn[nx*ny] = {0};
	for (int i = 0; i <nx; i++) {
		x[i] = i * 0.05;
		y[i] = i * 0.05;
	}

	  ofstream ufile("u.dat");
	  ofstream vfile("v.dat");
	  ofstream pfile("p.dat");

	
	for (int n = 0; n < nt; n++) {

		for (int j = ny; j < (nx-1)*(ny); j++){
			if (j % ny == 0 || j % ny == 40) {
				
			}
			else {
				b[j] = rho*((1/dt) * ((u[j+1] - u[j-1]) / (2*dx) + (v[j+nx] - v[j-nx]) / (2*dy)) - pow(((u[j+1]-u[j-1]) / (2*dx)), 2) -\
				2*(((u[j+nx] - u[j-nx]) / (2*dy)) * (v[j+1] - v[j-1]) / (2*dx)) - pow(((v[j+nx] - v[j-nx]) / (2*dy)),2));
			}
			if (j > (ny*(ny-3)) && j < (ny*(ny-2)))printf("b: %f uip: %f uim: %f ujp: %f ujm: %f vip: %f vim: %f vjp %f vjm %f id: %d |  ", b[j], u[j+1], u[j-1], u[j+nx], u[j-nx], v[j+1], v[j-1], v[j+nx], v[j-nx], j%ny);
		}
		


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
		for (int row = 0; row < (nx)*(ny); row++) {
			un[row] = u[row];
			vn[row] = v[row];	
		}
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

		printf("%34.30f \n", v[(20*ny)+3]);
		for(int row = 0; row < nx; row++) {
			u[row] = u[row * ny] = u[row * ny + (ny-1)] = 0.0;
			u[ny*(nx-1) + row] = 1.0;
			v[row] = v[ny*(nx-1)+row] = v[row*ny] = v[row*ny + (ny-1)] = 0.0;
		}
	
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
  ufile.close();
  vfile.close();
  pfile.close();



}
