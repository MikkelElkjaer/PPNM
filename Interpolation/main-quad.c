#include<stdio.h>
#include<stdlib.h>
#include<math.h>
typedef struct {int n; double *x, *y, *b, *c;} qspline;
qspline * qspline_alloc(int, double *, double*);
double qspline_eval(qspline *, double);
double qspline_derivative(qspline *, double);
double qspline_integral(qspline *, double);
void qspline_free(qspline *);

int main(int argc, char** argv) {
	/* make data file (should not be rerun): */
/*	for(int i=0; i<5; i++) printf("%i %i \n", i+1, 1);
	printf("\n \n");
	for(int i=0; i<5; i++) printf("%i %i \n", i+1, i+1);
	printf("\n \n");
        for(int i=0; i<5; i++) printf("%i %i \n", i+1, (i+1)*(i+1));
*/

int ntot = atoi(argv[1]) - 4; /* number of data points, and -4 to take two double lineskips into account */
int samples = 3; /* number of data samples stored in data_quad.txt */
int n = ntot/samples;
int b_analytical[3] = {0, 1, 0}, c_analytical[3] = {0, 0, 1}; /* b[2] will be appropriately altered in the for-loop below */
double deriv_analytical[3] = {0, 1, 2}; /* gets properly altered below */
double integ_analytical_factors[3] = {1, 1./2, 1./3}; /* the analytical integrals are calculated below */
double integ_analytical[3] = {0, 0, 0};
fprintf(stderr,"#b_i b_analytical,i c_i c_analytical,i \n");
for(int j=1;j<=samples;j++){ // run the same procedure for each of the three samples
	double x[n], y[n];
	for(int i=0; i<n;i++){
		scanf("%lg", &x[i]); /* alternatively: x+i */
		scanf("%lg", &y[i]);
	}
	qspline *s = qspline_alloc(n,x,y);

	for(int i=0;i<n-1;i++){
		if(j==3) b_analytical[j-1] +=2;
		fprintf(stderr,"%g %i %g %i\n", s->b[i], b_analytical[j-1], s->c[i], c_analytical[j-1]);
	}
	fprintf(stderr,"\n \n");

	for(double z=x[0];z<=x[n-1];z+=0.25){
		if(j==3){deriv_analytical[j-1]=2; deriv_analytical[j-1] *=z;} /* adjust analytical derivative */
		integ_analytical[j-1] = integ_analytical_factors[j-1]*pow(z,j) - integ_analytical_factors[j-1]*1; /* calculate the analytical integral */
		printf("%g %g %g %g %g %g\n",
			z, qspline_eval(s,z),
			qspline_derivative(s,z), deriv_analytical[j-1],
			qspline_integral(s,z), integ_analytical[j-1]);
	}
	qspline_free(s);
	printf("\n \n");
}
return 0;
}
