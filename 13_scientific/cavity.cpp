#include <cstdio>
#include <math.h>
#include <iostream>
using namespace std;
int main() {
	const int nx = 41;
       	const int ny = 41;
	const int nt = 50;
	double test = 0;
	const int nit = 50;
	const double dx = 2.0 / (nx - 1);
	const double dy = 2.0 / (ny-1);
	const double dt = 0.01;
	const int rho = 1;
	const double nu = 0.02;
	double x[nx];
	double y[ny];
	double u[nx][ny] = {0};
	double v[nx][ny] = {0};
	double p[nx][ny] = {0};
	double b[nx][ny] = {0};

	double pn[41][41],  un[41][41],  vn[41][41] = {};
	for (int i = 0; i <nx; i++) {
		x[i] = i * 0.05;
		y[i] = i * 0.05;
	}
	for (int n = 0; n < nt; n++) {
		for (int j = 1; j < (ny-1); j++) {
			for (int i = 1; i < (nx-1); i++) {
				b[j][i] = rho * ((1 / dt) *\
                    ((u[j][i+1] - u[j][i-1]) / (2 * dx) + (v[j+1][i] - v[j-1][i]) / (2 * dy)) - pow(((u[j][i+1] - u[j][i-1]) / (2 * dx)),2) - 2 * (((u[j+1][i] - u[j-1][i]) / (2 * dy)) *\
                     (v[j][i+1] - v[j][i-1]) / (2 * dx)) - pow(((v[j+1][i] - v[j-1][i]) / (2 * dy)),2));
	//		printf("%f \n", b[j][i]  ); 
		
//		printf("%f \n",(0-0)/(2*dx));	
			}
		}

	
		for (int it = 0; it < nit; it++) {
			for (int p1 = 0; p1 < nx; p1++) {
				for (int p2 = 0; p2 < ny; p2++) {
					pn[p1][p2] = p[p1][p2];
				}
			}
			for (int j = 1; j < (ny-1); j++) {
				for (int i = 1; i < (nx-1); i++) {
					p[j][i] = (pow(dy,2) * (pn[j][i+1] + pn[j][i-1]) +\
                          pow(dx,2) * (pn[j+1][i] + pn[j-1][i]) -  b[j][i] * pow(dx,2) * pow(dy,2)) / (2 * (pow(dx,2) + pow(dy,2)));
//	printf("%d \n", b[j][i]);
				}
			}
			for (int row = 0; row < nx; row++) {
				p[row][ny-1] = p[row][ny-2];
				p[0][row] = p[1][row];
				p[row][0] = p[row][1];
				p[nx-1][row] = 0;
			}
		}
		for (int row = 0; row < nx; row++) {
			for (int col = 0; col < ny; col++) {
				un[row][col] = u[row][col];
				vn[row][col] = v[row][col];
			}
		}
		for (int j = 1; j < ny-1; j++) {
			for (int i = 1; i < nx-1; i++) {
				u[j][i] = un[j][i] - un[j][i] * dt / dx * (un[j][i] - un[j][i - 1])\
                               - un[j][i] * dt / dy * (un[j][i] - un[j - 1][i])\
                               - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])\
                               + nu * dt / pow(dx,2) * (un[j][i+1] - 2 * un[j][i] + un[j][i-1])\
                               + nu * dt / pow(dy,2) * (un[j+1][i] - 2 * un[j][i] + un[j-1][i]);

				v[j][i] = vn[j][i] - vn[j][i] * dt / dx * (vn[j][i] - vn[j][i - 1])\
                               - vn[j][i] * dt / dy * (vn[j][i] - vn[j - 1][i])\
                               - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])\
                               + nu * dt / pow(dx,2) * (vn[j][i+1] - 2 * vn[j][i] + vn[j][i-1])\
                               + nu * dt / pow(dy,2) * (vn[j+1][i] - 2 * vn[j][i] + vn[j-1][i]);
			}
		}
		printf("%34.30f \n", v[20][2]);
	
		for (int row = 0; row < nx; row++) {
			u[0][row] = u[row][0] = u[row][ny-1] = 0;
			u[nx-1][row] = 1;
			v[0][row] = v[nx-1][row] = v[row][0] = v[row][ny-1] = 0;
		}
	
	}

}
