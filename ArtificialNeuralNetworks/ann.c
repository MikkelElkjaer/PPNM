#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_multimin.h>
#include"ann.h"
// data is kept in this form: data = {a_0,b_0,w_0,a_1,b_1,w_1,a_2,...}

ann* ann_alloc(int n, double (*f)(double)) {
        ann* network = malloc(sizeof(ann));
        network->n=n;
        network->f=f;
        network->data=gsl_vector_alloc(3*n);
        return network;
}

void ann_free(ann* network){
        gsl_vector_free(network->data);
        free(network);
}

double ann_feed_forward(ann* network, double x){
        double sum=0;
        for(int i=0;i < network->n;i++){
                double a = gsl_vector_get(network->data,3*i);
                double b = gsl_vector_get(network->data,3*i+1);
                double w = gsl_vector_get(network->data,3*i+2);
                sum += network->f((x-a)/b)*w;
        }
        return sum;
}

/* At the moment of solving this exercise I suspect there is a bug in my Quasi-Newton Minimizer,
and therefore, I will use the gsl-minimizer based on the Nelder-Mead routine instead. */

void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist){

	/* function to minimize: */
	double delta(const gsl_vector* p, void *params){
		gsl_vector_memcpy(network->data,p);
		double sum = 0;
		for(int i=0;i < xlist->size;i++){
			double x = gsl_vector_get(xlist,i);
			double ytable = gsl_vector_get(ylist,i);
			double y = ann_feed_forward(network,x);
			sum+=fabs(y-ytable);
		}
		return sum/xlist->size; // the factor of 1/xlist->size is mentioned at the lectures and is in principle optional
	}

	gsl_vector* p = gsl_vector_alloc(network->data->size);
	gsl_vector* h = gsl_vector_alloc(network->data->size);
	gsl_vector_memcpy(p,network->data); // p thus holds the parameters of the network. This acts as the initial guess
	for(int i=0;i < network->data->size;i++) gsl_vector_set(h,i, fabs(gsl_vector_get(p,i))/5); // set initial stepsizes of the simplex
	double size_goal = 1e-3;

	const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
	gsl_multimin_fminimizer *s = NULL;
	gsl_multimin_function minex_func;

	size_t iter = 0;
	int status; double size;
	double par[1] = {0}; // dummy params

	minex_func.n = network->data->size;
	minex_func.f = &delta;
	minex_func.params = (void *)par;

	s = gsl_multimin_fminimizer_alloc(T,network->data->size);
	gsl_multimin_fminimizer_set(s, &minex_func, p, h);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);
		if(status) break;

		size = gsl_multimin_fminimizer_size(s);
		status = gsl_multimin_test_size(size, size_goal);

	fprintf(stderr,"%li %g %g\n", iter, s->fval, size); // info about the minimization is printed to log.txt
	}
	while (status == GSL_CONTINUE && iter < 10000);


	gsl_vector_memcpy(network->data,s->x); // store the found parameters in the network
	gsl_vector_free(p); gsl_vector_free(h); gsl_multimin_fminimizer_free(s);
}
