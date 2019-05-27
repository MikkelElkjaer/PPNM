#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<math.h>

int ODE_driver(gsl_vector* tlist, gsl_matrix* ylist, double b, double h, double acc, double eps, int max,
		void stepper(double t, double h, gsl_vector* yt, void f(double t, gsl_vector* y, gsl_vector* dydt),
			     gsl_vector* yth, gsl_vector* err),
		void f(double t, gsl_vector* y, gsl_vector* dydt)) { // CHECK: type of y and dydt compared to exercise note...

	int j=0; const int n = ylist->size1;
	double t, err, normy, tol, a=gsl_vector_get(tlist,0);
	gsl_vector* y = gsl_vector_alloc(n);
	gsl_vector* yh = gsl_vector_alloc(n);
	gsl_vector* dy = gsl_vector_alloc(n);

	while(gsl_vector_get(tlist,j) < b) {
		t = gsl_vector_get(tlist,j);
		gsl_matrix_get_col(y,ylist,j);
		if(t+h>b) h = b-t; // avoid stepping too far

		stepper(t,h,y,f,yh,dy); // make a step and get yh and dy

		err = gsl_blas_dnrm2(dy);
		normy = gsl_blas_dnrm2(yh);

		tol = (normy*eps+acc)*sqrt(h/(b-a)); // as in eqn. (41) in the notes

		if(err<tol) { // that is, a valid step
			j++;
			if(j>max-1) return -1; // too many steps, return an indication hereof
			gsl_vector_set(tlist,j, t+h);
			gsl_matrix_set_col(ylist,j,yh);
		}
		if(err>0) h*=pow(tol/err,0.25)*0.95; else h*=2;
	}
return j; // return the index of the last element of tlist and last column of ylist
}
