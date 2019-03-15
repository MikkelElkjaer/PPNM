#include<gsl/gsl_multimin.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<gsl/gsl_vector.h>
#define TYPE gsl_multimin_fminimizer_nmsimplex2

struct expdata {int n; double *time, *activity, *error;};

double function_to_minimize (const gsl_vector *x, void *params) {
	double  A = gsl_vector_get(x,0);
	double  T = gsl_vector_get(x,1);
	double  B = gsl_vector_get(x,2);
	struct expdata *p = (struct expdata*) params;
	int     n = p->n;
	double *t = p->time;
	double *y = p->activity;
	double *e = p->error;
	double sum=0;
	#define f(t) A*exp(-(t)/T) + B
	for(int i=0;i<n;i++) sum += pow( (f(t[i]) - y[i])/e[i] ,2);
	return sum;
}

int main(int argc, char** argv){
	FILE* lsfstream = fopen("out_minimization.txt","w");
	int n=atoi(argv[1]);
	double activity[n],time[n],error[n];
	for(int i=0;i<n;i++)
		scanf("%lg %lg %lg",time+i,activity+i,error+i);

	int dim=3;
	struct expdata data;
	data.n=n;
	data.activity=activity;
	data.time=time;
	data.error=error;
	gsl_multimin_function F;
	F.f = function_to_minimize;
	F.n = dim;
	F.params=(void*)&data;

	gsl_multimin_fminimizer *M;
	M = gsl_multimin_fminimizer_alloc(TYPE,dim);
	gsl_vector* start=gsl_vector_alloc(dim);
	gsl_vector* step=gsl_vector_alloc(dim);
	gsl_vector_set(start,0,6); /* A, T, B */
	gsl_vector_set(start,1,2);
	gsl_vector_set(start,2,1);
	gsl_vector_set(step,0,0.5);
	gsl_vector_set(step,1,0.3);
	gsl_vector_set(step,2,0.3);

	gsl_multimin_fminimizer_set(M,&F,start,step);

	int iter=0,status;
	double size;
	do{
		iter++;
		status = gsl_multimin_fminimizer_iterate(M);
		if (status) break;

		size = gsl_multimin_fminimizer_size (M);
		status = gsl_multimin_test_size (size, 1e-2);

		if (status == GSL_SUCCESS)
          		fprintf (lsfstream,"converged to minimum at\n");

      	fprintf (lsfstream,"iter=%5i A=%8lg T=%8lg B=%8lg function_to_minimize=%8lg size=%8lg\n",
              iter,
              gsl_vector_get (M->x, 0),
              gsl_vector_get (M->x, 1),
              gsl_vector_get (M->x, 2),
              M->fval, size);

		if (status == GSL_SUCCESS)
		        fprintf (lsfstream,"\nThe estimated lifetime is thus: %lg \n", gsl_vector_get(M->x,1));

    	}
  while (status == GSL_CONTINUE && iter < 100);

	double A=gsl_vector_get(M->x,0);
	double T=gsl_vector_get(M->x,1);
	double B=gsl_vector_get(M->x,2);
	for(int i=0;i<n;i++)
		fprintf(stderr,"%g %g %g %g\n",
		*(time+i),*(activity+i),*(error+i)
		,A*exp(-time[i]/T)+B
		);

gsl_multimin_fminimizer_free(M);
gsl_vector_free(start);
gsl_vector_free(step);
fclose(lsfstream);
return 0;
}
