#include<gsl/gsl_vector.h>

void rkstep23(double t, double h, gsl_vector* yt, void f(double t, gsl_vector* y, gsl_vector* dydt),
		gsl_vector* yth, gsl_vector* err) { // ### y and yt the same perhaps?

	int i; const int n = yt->size;
	gsl_vector* k1 = gsl_vector_alloc(n);
	gsl_vector* k2 = gsl_vector_alloc(n);
	gsl_vector* k3 = gsl_vector_alloc(n);
	gsl_vector* k4 = gsl_vector_alloc(n);
	gsl_vector* yaux = gsl_vector_alloc(n);

	f(t,yt,k1); for(i=0;i<n;i++) gsl_vector_set(yaux,i, gsl_vector_get(yt,i) + 1./2*gsl_vector_get(k1,i)*h );
	f(t+1./2*h,yaux,k2); for(i=0;i<n;i++) gsl_vector_set(yaux,i, gsl_vector_get(yt,i) + 3./4*gsl_vector_get(k2,i)*h );
	f(t+3./4*h,yaux,k3); for(i=0;i<n;i++) gsl_vector_set(yth,i,
						gsl_vector_get(yt,i) + (2./9*gsl_vector_get(k1,i)
						+ 1./3*gsl_vector_get(k2,i) + 4./9*gsl_vector_get(k3,i))*h ); // this sets the first output: yth
	f(t+h,yth,k4);
	for(i=0;i<n;i++){
				gsl_vector_set(yaux,i, gsl_vector_get(yt,i) + (7./24*gsl_vector_get(k1,i)
					+ 1./4*gsl_vector_get(k2,i) + 1./3*gsl_vector_get(k3,i)
					+ 1./8*gsl_vector_get(k4,i))*h );
				gsl_vector_set(err,i, gsl_vector_get(yth,i) - gsl_vector_get(yaux,i)); // this sets the second output: err
			}


gsl_vector_free(k1); gsl_vector_free(k2); gsl_vector_free(k3); gsl_vector_free(k4); gsl_vector_free(yaux);
// HUSK FREE
}
