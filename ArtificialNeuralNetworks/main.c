#include<math.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>
#include"ann.h"

double activation_function(double x) { return x*exp(-x*x);}
double function_to_fit(double x) { return cos(4*x)*exp(-x*x);}

int main(){
	int n=10; // number of neurons
	ann* network = ann_alloc(n,&activation_function);

	/* make tabulated data: */
	double a=-M_PI/2, b = M_PI/2, x, y;
	int N = 30;
	gsl_vector* xlist = gsl_vector_alloc(N);
	gsl_vector* ylist = gsl_vector_alloc(N);
	for(int i=0;i<N;i++){
		x = a+(b-a)*i/(N-1);
		y = function_to_fit(x);
		gsl_vector_set(xlist,i,x);
		gsl_vector_set(ylist,i,y);
		printf("%g %g\n",x,y); // send table to output file
	}

	for(int i=0;i<n;i++){ // initialize parameters of the neurons
		gsl_vector_set(network->data,3*i, a+(b-a)*i/(N-1));
		gsl_vector_set(network->data,3*i+1, 1);
		gsl_vector_set(network->data,3*i+2, 1);
	}

	ann_train(network,xlist,ylist);

	printf("\n\n");

	double step = (b-a)/64;
	for(double z=a;z<=b;z+=step){
		y = ann_feed_forward(network,z);
		printf("%g %g\n",z,y);
	}

gsl_vector_free(xlist);gsl_vector_free(ylist);ann_free(network);
return 0;
}
