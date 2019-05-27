#include<math.h>
#include<assert.h>
#include<stdio.h>

double Integrator(double f(double), double, double, double, double, double*);

int main(){
/* nested functions will be exploited all throughout this main file */
/* test function 1: */
	int fcalls=0;
	double a=0, b=1, acc=1e-4, eps=1e-4; // limits of integration and absolute and relative accuracy demanded
	double err=0; // holds the error estimate of the integrator
	double f1(double x){
		fcalls++;
		return sqrt(x);};
	double Q = Integrator(f1,a,b,acc,eps,&err);
	printf("The integral of sqrt(x) on [0,1] is found to be: %g\nThe error estimate is: %g\n", Q,err);
	printf("The number of function calls used: %i\n", fcalls);

/* test function 2: */
	fcalls=0; err=0;
	double f2(double x){ fcalls++; return 1/sqrt(x);}
	Q = Integrator(f2,a,b,acc,eps,&err);
        printf("The integral of 1/sqrt(x) on [0,1] is found to be: %g\nThe error estimate is: %g\n", Q,err);
        printf("The number of function calls used: %i\n", fcalls);

/* test function 3: */
	fcalls=0; err=0;
	double f3(double x){ fcalls++; return log(x)/sqrt(x);}
	Q = Integrator(f3,a,b,acc,eps,&err);
        printf("The integral of ln(x)/sqrt(x) on [0,1] is found to be: %g\nThe error estimate is: %g\n", Q,err);
        printf("The number of function calls used: %i\n", fcalls);

/* test function 4: */
	fcalls=0; err=0;
	acc = 1e-13; eps = 1e-13; /* decreasing these to 1e-14 wont improve the estimate of pi
	(as we are limited to double precision) - rather, it only increases fcalls */
	double f4(double x){ fcalls++; return 4*sqrt(1-(1-x)*(1-x));}
	Q = Integrator(f4,a,b,acc,eps,&err);
        printf("The integral of 4*sqrt(1-(1-x)*(1-x)) on [0,1] is found to be: %.21g\nThe error estimate is: %g\n", Q,err);
        printf("The number of function calls used: %i\n", fcalls);

return 0;
}
