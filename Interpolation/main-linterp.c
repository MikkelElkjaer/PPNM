#include<stdio.h>
#include<stdlib.h>

double linterp(int, double *, double *, double);
double linterp_integ(int, double *, double *, double);

int main(int argc, char** argv){
	/* Linear Interpolation: */
	/* to make the tabulated data in data_linear.txt (should not be needed to run again): */
/*	for(double x=-10;x<=10;x+=0.5)
		printf("%g %g \n", x, 0.5*x*x); */

	int n = atoi(argv[1]); /* number of data points */
	double x[n], y[n];
	for(int i=0; i<n;i++){
		scanf("%lg", &x[i]); /* alternatively: x+i */
		scanf("%lg", &y[i]);
	}

	for(double z=x[0]; z<=x[n-1]; z+=0.2){
		double fz = linterp(n,x,y,z);
		double residual = 0.5*z*z - fz;
		double integ = linterp_integ(n,x,y,z);
		double integ_analytical = 1./6 * (z*z*z - x[0]*x[0]*x[0]);
		double integ_residual = integ_analytical - integ;
		/* quite naturally, both residuals depend heavily on the stepsize in the tabulated data... */
		printf("%g %g %g %g %g %g\n", z, fz, residual, integ, integ_analytical, integ_residual);
	}
return 0;
}
