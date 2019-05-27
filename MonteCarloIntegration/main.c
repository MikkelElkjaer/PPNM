#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

void plainmc(int, gsl_vector*, gsl_vector*, double f(gsl_vector*), int, double*, double*);

int main(){
/* test function 1 (a simple one to make sure it works as desired): */
	int dim = 2, N; double result, error;
	gsl_vector* a1 = gsl_vector_alloc(dim);
	gsl_vector* b1 = gsl_vector_alloc(dim);

	gsl_vector_set(a1,0,0); gsl_vector_set(b1,0,1);
	gsl_vector_set(a1,1,0); gsl_vector_set(b1,1,2);

	double f1(gsl_vector* x) { return gsl_vector_get(x,0)*gsl_vector_get(x,0) + sqrt(gsl_vector_get(x,1));}

	N = 10000;
	plainmc(dim,a1,b1,f1,N,&result,&error);

	printf("Using %i points we get for test function 1:\nResult = %g\nError = %g\n", N, result, error);
	printf("The actual result is 2.55228\n");

/* test function 2 (1s Hydrogen wavefunction within 5 Bohr radii (a_0=1)):
NOTE that the boundaries of the integration are chosen to cover only the 1/8th part of real space
where x,y,z are all positive. We simply multiply the result by 8 (or equivalently the integrand)
to get the right result (something close to 1) as we exploit the spherical symmetry of the 1s wavefct.
We perform this trick to achieve a better resolution of the integrand.*/
	dim = 3;
	gsl_vector* a2 = gsl_vector_alloc(dim);
	gsl_vector* b2 = gsl_vector_alloc(dim);

	gsl_vector_set(a2,0,0); gsl_vector_set(b2,0,5); // boundaries
	gsl_vector_set(a2,1,0); gsl_vector_set(b2,1,5);
	gsl_vector_set(a2,2,0); gsl_vector_set(b2,2,5);

	double f2(gsl_vector* x) { // the Bohr radius a_0 is set to 1
		double r = sqrt(gsl_vector_get(x,0)*gsl_vector_get(x,0)
			+ gsl_vector_get(x,1)*gsl_vector_get(x,1) + gsl_vector_get(x,2)*gsl_vector_get(x,2));
		return 1/M_PI*exp(-2*r)*8;
	}
	N = 1000000;
	plainmc(dim,a2,b2,f2,N,&result,&error);

        printf("Using %i points we get for test function 2:\nResult = %g\nError = %g\n", N, result, error);
        printf("The result should be close to (and less than) 1\n");

/* test function 3 (the difficult singular one from the exercise): */
	dim = 3;
	gsl_vector* a3 = gsl_vector_alloc(dim); gsl_vector* b3 = gsl_vector_alloc(dim);

	gsl_vector_set(a3,0,0); gsl_vector_set(b3,0,M_PI); // boundaries
        gsl_vector_set(a3,1,0); gsl_vector_set(b3,1,M_PI);
        gsl_vector_set(a3,2,0); gsl_vector_set(b3,2,M_PI);

	double f3(gsl_vector* x) {
		return 1/(M_PI*M_PI*M_PI)*1/(1 - cos(gsl_vector_get(x,0))*cos(gsl_vector_get(x,1))*cos(gsl_vector_get(x,2)));
	}
	N = 1000000;
	plainmc(dim,a3,b3,f3,N,&result,&error);
        printf("Using %i points we get for the highly singular test function 3:\nResult = %g\nError = %g\n", N, result, error);
        printf("The result should be approx. 1.39320393\n");

/* PART B - error behaves as O(1/sqrt(N)): */
/* We will reuse test function 1 for simplicity */
	dim = 2; int iters = 10, Nmax = 1000000; // (1 milllion)
	double errorsum;

for(N=100;N<Nmax+1;N*=10){
	errorsum = 0;
	for(int i=0;i<iters;i++){
		plainmc(dim,a1,b1,f1,N,&result,&error);
		errorsum += error;
	}
	errorsum/=iters;
	fprintf(stderr,"%i %g\n",N,errorsum);
}

/* the resulting plot shows decent agreement with the expected 1/sqrt(N)-scaling */

return 0;
}
