#include<assert.h>

double linterp_integ(int n, double *x, double *y, double z){
	assert(n>1 && z>=x[0] && z<=x[n-1]); /* make sure that z is in fact in spline-able range */
        /* do binary search: */
        int lower = 0, upper = n-1;
        while(upper-lower>1){
                int mid = (upper+lower)/2;
                if(z>x[mid]) lower=mid;
                else upper=mid;
        }
	double integ = 0;
	for(int i=0;i<lower;i++)
		integ += (x[i+1]-x[i]) * (y[i] + 0.5*(y[i+1]-y[i]));

	double fz = (y[upper]-y[lower])/(x[upper]-x[lower]) * (z-x[lower]) + y[lower];
	integ += (z-x[lower]) * (y[lower] + 0.5*(fz-y[lower]));
return integ;
}
