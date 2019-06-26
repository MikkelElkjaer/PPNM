#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>
#include<gsl/gsl_complex.h>
#include<gsl/gsl_complex_math.h>

int ODE_driver(gsl_vector_complex* tlist, gsl_matrix_complex* ylist, gsl_complex b, gsl_complex h, double acc, double eps, int max,
		void stepper(gsl_complex t, gsl_complex h, gsl_vector_complex* yt,
			void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt),
			gsl_vector_complex* yth, gsl_vector_complex* err),
		void f(gsl_complex t, gsl_vector_complex* y, gsl_vector_complex* dydt)) {

	int j=0; const int n = ylist->size1;
	gsl_complex t, a = gsl_vector_complex_get(tlist,0);
	double err, normy, tol, accumulatedstepsdist;
	gsl_vector_complex* y = gsl_vector_complex_alloc(n);
	gsl_vector_complex* yh = gsl_vector_complex_alloc(n);
	gsl_vector_complex* dy = gsl_vector_complex_alloc(n);

	/* enforce that a straight line from a to b is taken: */
	double abdistance = gsl_complex_abs(gsl_complex_sub(b,a));
	double relength = GSL_REAL(b) - GSL_REAL(a);
	double imlength = GSL_IMAG(b) - GSL_IMAG(a);
	double realresolution = GSL_REAL(h)/relength, imagresolution = GSL_REAL(h)/imlength, smallestres;
	if(realresolution < imagresolution) smallestres = realresolution; else smallestres = imagresolution;
	h = gsl_complex_rect(smallestres*relength, smallestres*imlength);
	accumulatedstepsdist = 0;

	while(accumulatedstepsdist < (abdistance-0.0001*acc)) {
// the 0.0001*acc-factor is a hotfix to mitigate a somewhat rare event where the loop was not broken with the last step as the accumulatedstepdistance
// would still seem less than abdistance due to round off error at the order of the machine epsilon... I have not found a more elegant fix to this.
		t = gsl_vector_complex_get(tlist,j);
		gsl_matrix_complex_get_col(y,ylist,j);
		if((accumulatedstepsdist + gsl_complex_abs(h)) > abdistance) h = gsl_complex_sub(b,t); // avoid stepping too far

		stepper(t,h,y,f,yh,dy); // make a step and get yh and dy

		err = gsl_blas_dznrm2(dy);
		normy = gsl_blas_dznrm2(yh);

		tol = (normy*eps+acc)*sqrt(gsl_complex_abs(h)/abdistance); // as in eqn. (41) in the notes

		if(err<tol) { // that is, a valid step
			j++;
			accumulatedstepsdist += gsl_complex_abs(h);
			if(j>max-1) return -1; // too many steps, return an indication hereof
			gsl_vector_complex_set(tlist,j, gsl_complex_add(t,h));
			gsl_matrix_complex_set_col(ylist,j,yh);
		}
		h = gsl_complex_mul_real(h, pow(tol/err,0.25)*0.95); // update step as in eqn. (40)

	}
gsl_vector_complex_free(y); gsl_vector_complex_free(yh); gsl_vector_complex_free(dy);
return j; // return the index of the last element of tlist and last column of ylist
}
