#include<gsl/gsl_vector.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

void rkstep23(gsl_complex t, gsl_complex h, gsl_vector_complex* yt, void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt),
		gsl_vector_complex* yth, gsl_vector_complex* err) {

	int i; const int n = yt->size;
	double stepsize = gsl_complex_abs(h);
	gsl_complex z1, z2, z3, z4; // dummies
	gsl_vector_complex* k1 = gsl_vector_complex_alloc(n);
	gsl_vector_complex* k2 = gsl_vector_complex_alloc(n);
	gsl_vector_complex* k3 = gsl_vector_complex_alloc(n);
	gsl_vector_complex* k4 = gsl_vector_complex_alloc(n);
	gsl_vector_complex* yaux = gsl_vector_complex_alloc(n);

	f(t,yt,k1);
	for(i=0;i<n;i++) {
		z1 = gsl_complex_mul_real(gsl_vector_complex_get(k1,i),stepsize*1./2);
		z2 = gsl_complex_add(gsl_vector_complex_get(yt,i),z1);
		gsl_vector_complex_set(yaux,i,z2);
	}
	f(gsl_complex_add(t,gsl_complex_mul_real(h,1./2)),yaux,k2); //f(t+1./2*h,yaux,k2);
	for(i=0;i<n;i++) {
		z1 = gsl_complex_mul_real(gsl_vector_complex_get(k2,i),stepsize*3./4);
		z2 = gsl_complex_add(gsl_vector_complex_get(yt,i),z1);
		gsl_vector_complex_set(yaux,i,z2);
	}
	f(gsl_complex_add(t,gsl_complex_mul_real(h,3./4)),yaux,k3); //f(t+3./4*h,yaux,k3);
	for(i=0;i<n;i++) {
		z1 = gsl_complex_mul_real(gsl_vector_complex_get(k3,i),4./9);
		z2 = gsl_complex_mul_real(gsl_vector_complex_get(k2,i),1./3);
		z3 = gsl_complex_mul_real(gsl_vector_complex_get(k1,i),2./9);
		z2 = gsl_complex_add(z2,z3);
		z1 = gsl_complex_add(z1,z2);
		z3 = gsl_complex_mul_real(z1,stepsize);
		z2 = gsl_complex_add(gsl_vector_complex_get(yt,i),z3);
		gsl_vector_complex_set(yth,i,z2); // this sets the first output: yth
		}
	f(gsl_complex_add(t,h),yth,k4); //f(t+h,yth,k4);
	for(i=0;i<n;i++) {
		z1 = gsl_complex_mul_real(gsl_vector_complex_get(k4,i),1./8);
		z2 = gsl_complex_mul_real(gsl_vector_complex_get(k3,i),1./3);
		z3 = gsl_complex_mul_real(gsl_vector_complex_get(k2,i),1./4);
		z4 = gsl_complex_mul_real(gsl_vector_complex_get(k1,i),7./24);
		z3 = gsl_complex_add(z3,z4); z1 = gsl_complex_add(z1,z2);
		z1 = gsl_complex_add(z1,z3);
		z1 = gsl_complex_mul_real(z1,stepsize);
		z2 = gsl_vector_complex_get(yt,i);
		z3 = gsl_complex_add(z1,z2);
		gsl_vector_complex_set(yaux,i,z3);

		z1 = gsl_vector_complex_get(yth,i);
		z4 = gsl_complex_sub(z1,z3);
		gsl_vector_complex_set(err,i,z4); // this sets the second output: err
	}

gsl_vector_complex_free(k1); gsl_vector_complex_free(k2); gsl_vector_complex_free(k3); gsl_vector_complex_free(k4); gsl_vector_complex_free(yaux);
}
