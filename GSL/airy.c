#include<stdio.h>
#include<gsl/gsl_sf_airy.h>

void airy(void){
	for(double x=-15.0;x<5;x+=0.1){
	double a = gsl_sf_airy_Ai(x,GSL_PREC_DOUBLE);
	double b = gsl_sf_airy_Bi(x,GSL_PREC_DOUBLE);
	fprintf(stderr, "%g %g %g \n", x,a,b);
	}
}

