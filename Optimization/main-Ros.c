#include<gsl/gsl_multimin.h>
#include<stdio.h>

double Rosenbrock(const gsl_vector *vec, void *params){
	double x=gsl_vector_get(vec,0);
	double y=gsl_vector_get(vec,1);
	double *p = (double *)params;

	return (p[0]-x)*(p[0]-x) + p[1]*(y-x*x)*(y-x*x);
}

int main(void) {
	double par[2] = {1.0, 100.0};

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_vector *ss, *x;
	gsl_multimin_function minim_f;

	size_t iter = 0;
	int status;
	double size;

	/* Starting point */
	x = gsl_vector_alloc (2);
	gsl_vector_set (x, 0, 5);
	gsl_vector_set (x, 1, 5);

	/* Set initial step sizes to 0.1 */
	ss = gsl_vector_alloc (2);
	gsl_vector_set_all (ss, 0.1);

	minim_f.n = 2;
	minim_f.f = Rosenbrock;
	minim_f.params = par;

	s = gsl_multimin_fminimizer_alloc (T, 2);
	gsl_multimin_fminimizer_set (s, &minim_f, x, ss);

        printf("# iter, x, y, f(x,y), size\n");

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

      		if (status)
        		break;

      		size = gsl_multimin_fminimizer_size (s);
      		status = gsl_multimin_test_size (size, 1e-4);

      		if (status == GSL_SUCCESS)
        	{
          		printf ("converged to minimum at\n");
        	}

      		printf ("%5li %10.3lg %10.3lg %10.3lg %10.3lg\n",
              	iter,
              	gsl_vector_get (s->x, 0),
              	gsl_vector_get (s->x, 1),
              	s->fval, size);
    	}
	while (status == GSL_CONTINUE && iter < 500);

	gsl_vector_free(x);
	gsl_vector_free(ss);
	gsl_multimin_fminimizer_free (s);

	return status;
}
