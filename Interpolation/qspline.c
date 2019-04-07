#include<assert.h>
#include<stdlib.h>

typedef struct {int n; double *x, *y, *b, *c;} qspline;

qspline * qspline_alloc(int n, double *x, double *y){ // builds the qspline
	/* allocate: */
	qspline *s = (qspline*)malloc(sizeof(qspline));
	s->n = n;
	s->x = (double*)malloc(n*sizeof(double));
	s->y = (double*)malloc(n*sizeof(double));
	s->b = (double*)malloc((n-1)*sizeof(double));
	s->c = (double*)malloc((n-1)*sizeof(double));
	/* fill in with data: */
	for(int i=0;i<n;i++){
		s->x[i] = x[i];
		s->y[i] = y[i];
	}
	/* calculate spline parameters b_i and c_i: */
	double Deltax[n-1], p[n-1];
	for(int i=0;i<n-1;i++){Deltax[i] = x[i+1]-x[i]; p[i] = (y[i+1]-y[i])/Deltax[i];}
	s->c[0] = 0; // starting point of the upwards recursion via (11) in the notes
	for(int i=0;i<n-2;i++)
		s->c[i+1] = (p[i+1]-p[i]-s->c[i]*Deltax[i])/Deltax[i+1];
	s->c[n-2]/=2; // to setup downwards recursion
	for(int i=n-3;i>=0;i--)
		s->c[i] = (p[i+1]-p[i]-s->c[i+1]*Deltax[i+1])/Deltax[i];
	for(int i=0;i<n-1;i++)
		s->b[i] = p[i]-s->c[i]*Deltax[i];
return s; //note that the returned form is the one suitable for differentiation and integration
}

int binary_search(int n, double *x, double z){
        int lower = 0, upper = n-1, mid;
        while(upper-lower>1){
                mid = (upper+lower)/2;
                if(z>x[mid]) lower=mid;
                else upper=mid;
        }
return lower;
}

double qspline_eval(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]); // defensive assertion: make sure z is in splined x-range
	int i = binary_search(s->n, s->x,z);
	double Deltax = z - s->x[i];
return s->y[i] + Deltax*s->b[i] + Deltax*Deltax*s->c[i];
}

void qspline_free(qspline *s){
	free(s->x); free(s->y); free(s->b); free(s->c); free(s);
}

double qspline_derivative(qspline *s, double z){
	assert(z>=s->x[0] && z<=s->x[s->n-1]);
	int i = binary_search(s->n, s->x,z);
return s->b[i] + 2*s->c[i]*(z-s->x[i]);
}

double qspline_integral(qspline *s, double z){
        assert(z>=s->x[0] && z<=s->x[s->n-1]);
        int i = binary_search(s->n, s->x, z);
	double integral = 0;
	double dx;
	if(i>0){
		for(int j=0; j<i; j++){
			dx = s->x[j+1] - s->x[j];
			integral += s->y[j]*dx + s->b[j]*dx*dx/2 + s->c[j]*dx*dx*dx/3;
		}
	}
	dx = z - s->x[i];
	assert(dx<= s->x[i+1] - s->x[i]); // assert that binary_search found the right interval
	integral += s->y[i]*dx + s->b[i]*dx*dx/2 + s->c[i]*dx*dx*dx/3;
return integral;
}
