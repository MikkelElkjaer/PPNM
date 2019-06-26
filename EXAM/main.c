#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

int ODE_driver(gsl_vector_complex*, gsl_matrix_complex*, gsl_complex, gsl_complex, double, double, int,
		void stepper(gsl_complex, gsl_complex, gsl_vector_complex*, void f(gsl_complex, gsl_vector_complex*, gsl_vector_complex*),
				gsl_vector_complex*, gsl_vector_complex*), void f(gsl_complex, gsl_vector_complex*, gsl_vector_complex*));
void rkstep23(gsl_complex, gsl_complex, gsl_vector_complex*, void f(gsl_complex, gsl_vector_complex*, gsl_vector_complex*),
		 gsl_vector_complex*, gsl_vector_complex*);

void ode_real_exp(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt){
	gsl_vector_complex_set(dydt,0,gsl_vector_complex_get(y,0));
}

void ode_imag_exp(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt){
	gsl_vector_complex_set(dydt,0, gsl_complex_mul_real(gsl_vector_complex_get(y,0),-1));
}

void ode_complex_exp(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt){
	gsl_vector_complex_set(dydt,0, gsl_complex_mul_imag(gsl_vector_complex_get(y,0),1));
}

int main() {
/* test 1. Integrating exp(x) of a real variable. This is done to check that the generalization
to the complex plane did not accidentally introduce errors in the real domain. */
int max = 50; // maximum number of allowed steps for the driver to take

gsl_vector_complex* tlist = gsl_vector_complex_alloc(max);
gsl_matrix_complex* ylist = gsl_matrix_complex_alloc(1,max);
gsl_complex a, b, hstart, yinit = gsl_complex_rect(1,0);
double acc = 1e-4, eps = 1e-4;
GSL_SET_COMPLEX(&a,0,0); GSL_SET_COMPLEX(&b,3,0);
GSL_SET_COMPLEX(&hstart,0.1,0.1); // the driver should recognize that this hstart won't result in a straight line and correct it!
gsl_vector_complex_set(tlist,0,a); gsl_matrix_complex_set(ylist,0,0,yinit); // initial conditions exp(0+0*i)=1+0*i

int idxmax = ODE_driver(tlist, ylist, b, hstart, acc, eps, max, rkstep23, ode_real_exp); // apply driver. 1 means success, -1 means failure
if(idxmax<0) fprintf(stderr,"Driver failed in test 1. Try increasing the number of allowed steps.");


for(int i=0;i<idxmax+1;i++) fprintf(stderr,"%7.4g %7.4g %7.4g %7.4g %7.4g\n",
					GSL_REAL(gsl_vector_complex_get(tlist,i)), GSL_IMAG(gsl_vector_complex_get(tlist,i)),
					GSL_REAL(gsl_matrix_complex_get(ylist,0,i)), GSL_IMAG(gsl_matrix_complex_get(ylist,0,i)),
					exp(GSL_REAL(gsl_vector_complex_get(tlist,i))));

fprintf(stderr,"\n\n");
/* test 2. Integrating the complex exponential along the imaginary axis: */
max = 250;

gsl_vector_complex* tlist2 = gsl_vector_complex_alloc(max);
gsl_matrix_complex* ylist2 = gsl_matrix_complex_alloc(1,max);

acc = 1e-5; eps = 1e-5;
GSL_SET_COMPLEX(&a,0,0); GSL_SET_COMPLEX(&b,0,2*M_PI);
yinit = gsl_complex_rect(1,0);
GSL_SET_COMPLEX(&hstart,0.1,0.1); // the driver should again correct this poor guess
gsl_vector_complex_set(tlist2,0,a); gsl_matrix_complex_set(ylist2,0,0,yinit);

idxmax = ODE_driver(tlist2, ylist2, b, hstart, acc, eps, max, rkstep23, ode_imag_exp);
if(idxmax<0) fprintf(stderr,"Driver failed in test 2. Try increasing the number of allowed steps.");

for(int i=0;i<idxmax+1;i++) fprintf(stderr,"%7.5g %7.5g %7.5g %7.5g %7.5g\n",
					GSL_REAL(gsl_vector_complex_get(tlist2,i)), GSL_IMAG(gsl_vector_complex_get(tlist2,i)),
                                        GSL_REAL(gsl_matrix_complex_get(ylist2,0,i)), GSL_IMAG(gsl_matrix_complex_get(ylist2,0,i)),
					exp(-1*GSL_IMAG(gsl_vector_complex_get(tlist2,i))));

fprintf(stderr,"\n\n");

/* test 3: Integrating the complex exponential along the real axis: */
max = 300;
gsl_vector_complex* tlist3 = gsl_vector_complex_alloc(max);
gsl_matrix_complex* ylist3 = gsl_matrix_complex_alloc(1,max);
GSL_SET_COMPLEX(&a,0,0); GSL_SET_COMPLEX(&b,2*M_PI,0);
yinit = gsl_complex_rect(1,0);
GSL_SET_COMPLEX(&hstart,0.1,0.1);
gsl_vector_complex_set(tlist3,0,a); gsl_matrix_complex_set(ylist3,0,0,yinit);

idxmax = ODE_driver(tlist3, ylist3, b, hstart, acc, eps, max, rkstep23, ode_complex_exp);
if(idxmax<0) fprintf(stderr,"Driver failed in test 3. Try increasing the number of allowed steps.");

for(int i=0;i<idxmax+1;i++) fprintf(stderr,"%7.5g %7.5g %7.5g %7.5g %7.5g %7.5g\n",
                                        GSL_REAL(gsl_vector_complex_get(tlist3,i)), GSL_IMAG(gsl_vector_complex_get(tlist3,i)),
                                        GSL_REAL(gsl_matrix_complex_get(ylist3,0,i)), GSL_IMAG(gsl_matrix_complex_get(ylist3,0,i)),
					cos(GSL_REAL(gsl_vector_complex_get(tlist3,i))), sin(GSL_REAL(gsl_vector_complex_get(tlist3,i))));


gsl_vector_complex_free(tlist); gsl_vector_complex_free(tlist2); gsl_matrix_complex_free(ylist); gsl_matrix_complex_free(ylist2);
gsl_vector_complex_free(tlist3); gsl_matrix_complex_free(ylist3);
return 0;
}
