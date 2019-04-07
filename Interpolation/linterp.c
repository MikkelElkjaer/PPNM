#include<assert.h>

double linterp(int n, double *x, double *y, double z){
	assert(n>1 && z>=x[0] && z<=x[n-1]); /* make sure that z is in fact in spline-able range */
	/* do binary search: */
	int lower = 0, upper = n-1, mid;
	while(upper-lower>1){
		mid = (upper+lower)/2;
		if(z>x[mid]) lower=mid;
		else upper=mid;
	}
	return (y[upper]-y[lower])/(x[upper]-x[lower]) * (z-x[lower]) + y[lower]; /* a*(z-x_left) + b */
}

