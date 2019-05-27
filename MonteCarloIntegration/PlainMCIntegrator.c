#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_vector.h>
#include<time.h>
#define RND ((double)rand()/RAND_MAX)

void randompoint(int dim, gsl_vector* a, gsl_vector* b, gsl_vector* x) {
	for(int i=0;i<dim;i++) gsl_vector_set(x,i, gsl_vector_get(a,i)+RND*(gsl_vector_get(b,i)-gsl_vector_get(a,i)));
}

void plainmc(int dim, gsl_vector* a, gsl_vector* b, double f(gsl_vector* x), int N, double* res, double* err) {
	srand(time(0));
	double V=1;
	for(int i=0;i<dim;i++) V*=gsl_vector_get(b,i) - gsl_vector_get(a,i); // make proper volume
	double fsum=0, fsquaresum=0, fx;
	gsl_vector* x = gsl_vector_alloc(dim);

	for(int i=0;i<N;i++) {
		randompoint(dim,a,b,x);
		fx = f(x);
		fsum += fx;
		fsquaresum += fx*fx;
	}
	double avg = fsum/N;
	double var = fsquaresum/N - avg*avg;

	*res = avg*V;
	*err = sqrt(var/N)*V;

gsl_vector_free(x);
}
