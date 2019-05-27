#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>

int ODE_driver(gsl_vector*, gsl_matrix*, double, double, double, double, int,
		void stepper(double, double, gsl_vector*, void f(double, gsl_vector*, gsl_vector*),
				gsl_vector*, gsl_vector*), void f(double, gsl_vector*, gsl_vector*));
void rkstep23(double, double, gsl_vector*, void f(double, gsl_vector*, gsl_vector*), gsl_vector*, gsl_vector*);


void ode_cubic_root(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0, 1./(3*gsl_vector_get(y,0)*gsl_vector_get(y,0)));
}

void ode_Lidocaine(double t, gsl_vector* y, gsl_vector* dydt){
	gsl_vector_set(dydt,0, -0.09*gsl_vector_get(y,0) + 0.038*gsl_vector_get(y,1)); // dx/dt
	gsl_vector_set(dydt,1, 0.066*gsl_vector_get(y,0) - 0.038*gsl_vector_get(y,1)); // dy/dt
}

int main() {
/* first (very simple) test is the ode-representation of the cubic root: */
int max = 50; // maximum number of allowed steps for the driver to take
gsl_vector* tlist = gsl_vector_alloc(max);
gsl_matrix* ylist = gsl_matrix_alloc(1,max);
double a=1, b=10, hstart=0.1, acc=1e-5, eps=1e-5;

gsl_vector_set(tlist,0,a); gsl_matrix_set(ylist,0,0,1); // initial conditions y(1)=1

int idxmax = ODE_driver(tlist, ylist, b, hstart, acc, eps, max, rkstep23, ode_cubic_root); // apply driver. 1 means success, -1 means failure
if(idxmax<0) fprintf(stderr,"Driver failed. Try increasing the number of allowed steps.");

for(int i=0;i<idxmax+1;i++) fprintf(stderr,"%g %g\n", gsl_vector_get(tlist,i), gsl_matrix_get(ylist,0,i));

/* Second test: Lidocaine for patients with irregular heartbeats.
See details and short description at http://www.math.utah.edu/~gustafso/2250systems-de.pdf (visited 22/05/19) */

max = 50;
gsl_vector* tlist2 = gsl_vector_alloc(max);
gsl_matrix* ylist2 = gsl_matrix_alloc(2,max);

a = 0; b = 10; hstart = 0.1; acc = 1e-5; eps = 1e-5;

double initial_dosage = 2.0;
gsl_vector_set(tlist2,0, 0);
gsl_matrix_set(ylist2,0,0, 0); gsl_matrix_set(ylist2,1,0, initial_dosage); // initial conditions

idxmax = ODE_driver(tlist2,ylist2,b,hstart,acc,eps,max,rkstep23,ode_Lidocaine);
if(idxmax<0) printf("Driver failed. Try increasing the number of allowed steps.");

for(int i=0;i<=idxmax;i++) fprintf(stdout,"%g %g %g\n", gsl_vector_get(tlist2,i), gsl_matrix_get(ylist2,0,i), gsl_matrix_get(ylist2,1,i));

return 0;
}
