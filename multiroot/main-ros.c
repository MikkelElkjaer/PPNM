#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multiroots.h>

struct rparams {double a; double b;};

void print_state (size_t iter, gsl_multiroot_fsolver * s, struct rparams * params) {
	double x = gsl_vector_get (s->x,0);
	double y = gsl_vector_get (s->x,1);
	double a = params->a;
	double b = params->b;
	printf ("%3li %5.3lg %5.3lg %5.3lg %6.3lg %7.3lg \n",
	iter, x, y,
        gsl_vector_get (s->f, 0),
        gsl_vector_get (s->f, 1),
	(a-x)*(a-x)+b*(y-x*x)*(y-x*x));

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

  return GSL_SUCCESS;
}

int main(void) {
	const gsl_multiroot_fsolver_type *T;
  	gsl_multiroot_fsolver *s;

  	int status;
	size_t iter = 0;

	const size_t n = 2;
  	struct rparams p = {1.0, 100.0};
  	gsl_multiroot_function f = {&rosenbrock_deriv, n, &p};

  	double x_init[2] = {1.5, 2.5};
	gsl_vector *x = gsl_vector_alloc (n);

  	gsl_vector_set (x, 0, x_init[0]);
  	gsl_vector_set (x, 1, x_init[1]);

  	T = gsl_multiroot_fsolver_hybrids;
  	s = gsl_multiroot_fsolver_alloc (T, 2);
  	gsl_multiroot_fsolver_set (s, &f, x);

	printf("iter,  x,    y,  grad_x, grad_y, function_value \n");
	print_state (iter, s, &p);

  	do {
      		iter++;
      		status = gsl_multiroot_fsolver_iterate (s);

      		print_state (iter, s, &p);

      		if (status)   /* stop the solver if stuck */
        		break;

      		status = gsl_multiroot_test_residual (s->f, 1e-7);
    	}
  	while (status == GSL_CONTINUE && iter < 1000);

  	printf ("status = %s\n", gsl_strerror (status));

  	gsl_multiroot_fsolver_free (s);
  	gsl_vector_free (x);
  	return 0;
}

