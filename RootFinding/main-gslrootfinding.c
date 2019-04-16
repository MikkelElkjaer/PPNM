#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>

struct rparams {double a; double b; int fcalls;};
FILE *appendstream;

void print_state (size_t iter, gsl_multiroot_fsolver * s, struct rparams * params) {
	double x = gsl_vector_get (s->x,0);
	double y = gsl_vector_get (s->x,1);
	fprintf(appendstream,"iterations=%li, functions calls=%i, x=%.3lg, y=%.3lg \n",
	iter, params->fcalls, x, y);
}

int rosenbrock_deriv (const gsl_vector * x, void *params, gsl_vector * f) {
  double a = ((struct rparams *) params)->a;
  double b = ((struct rparams *) params)->b;
  const double x0 = gsl_vector_get (x, 0);
  const double x1 = gsl_vector_get (x, 1);
  const double y0 = -2*(a-x0)-4*b*x0*(x1-x0*x0); /* analytically calculated gradients */
  const double y1 = 2*b*(x1-x0*x0);

  gsl_vector_set (f, 0, y0);
  gsl_vector_set (f, 1, y1);
  ((struct rparams *) params)->fcalls++; // keep track of how many times the function is called
  return GSL_SUCCESS;}

int Himmelblau(const gsl_vector *x, void *params, gsl_vector *f){
        double a = 11, b = 7, y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(f,0, 2*(z*z+y-a)*2*z + 2*(z+y*y-b));
        gsl_vector_set(f,1, 2*(z*z+y-a) + 2*(z+y*y-b)*2*y);
	((struct rparams *) params)->fcalls++;
return GSL_SUCCESS;}

int main(void) {
	appendstream = fopen("out.txt","a");
	const gsl_multiroot_fsolver_type *T;
  	gsl_multiroot_fsolver *s;
  	int status;
	size_t iter = 0;
	const size_t n = 2;
  	struct rparams p = {1.0, 100.0, 0};
  	gsl_multiroot_function f = {&rosenbrock_deriv, n, &p};
  	double x_init[2] = {0.5, 0.5};
	gsl_vector *x = gsl_vector_alloc (n);
  	gsl_vector_set (x, 0, x_init[0]);
  	gsl_vector_set (x, 1, x_init[1]);

  	T = gsl_multiroot_fsolver_hybrids;
  	s = gsl_multiroot_fsolver_alloc (T, 2);
  	gsl_multiroot_fsolver_set (s, &f, x);

	fprintf(appendstream,"\nFrom the GSL-library routines:\n");
	fprintf(appendstream,"Rosenbrock minimum with the same intial guess and precision demands:\n");
  	do {
      		iter++;
      		status = gsl_multiroot_fsolver_iterate (s);
      		if (status)   /* stop the solver if stuck */
        		break;
      		status = gsl_multiroot_test_residual (s->f, 1e-5);
    	}
  	while (status == GSL_CONTINUE && iter < 1000);
	print_state(iter,s,&p); // for this exercise we only want the end-state

	/* Himmelblau part: */
	fprintf(appendstream,"\nHimmelblau minimum with the same initial guess and precision demands:\n");
	struct rparams q = {0,0,0};
	iter = 0;
	gsl_multiroot_fsolver *s2;
	gsl_multiroot_function g = {&Himmelblau, n, &q};
	gsl_vector_set(x,0,4.4); gsl_vector_set(x,1,4.2); // initial guess
	s2 = gsl_multiroot_fsolver_alloc (T, 2);
	gsl_multiroot_fsolver_set (s2, &g, x);

	do {
                iter++;
                status = gsl_multiroot_fsolver_iterate (s2);
                if (status)   /* stop the solver if stuck */
                        break;
                status = gsl_multiroot_test_residual (s2->f, 1e-5);
        }
        while (status == GSL_CONTINUE && iter < 1000);
        print_state(iter,s2,&q); // for this exercise we only want the end-state

  	gsl_multiroot_fsolver_free (s);
	gsl_multiroot_fsolver_free (s2);
  	gsl_vector_free (x);
  	return 0;
}

