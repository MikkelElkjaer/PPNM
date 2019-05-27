#include<math.h>
#include<assert.h>

double Integrator24(double f(double), double a, double b, double acc, double eps,
		    double f2, double f3, int recdepth, double* err) {

	assert(recdepth<1000000);
	double f1 = f(a+(b-a)/6);
	double f4 = f(a+5*(b-a)/6);
	double Q = (2*f1+f2+f3+2*f4)/6*(b-a); // higher order estimate
	double q = (f1+f2+f3+f4)/4*(b-a); // lower order estimate
	double tol = acc+eps*fabs(Q);
	double localerr = fabs(Q-q);

	if(localerr < tol){ *err += localerr; return Q;}
	else {
		double Q1=Integrator24(f,a,(a+b)/2,acc/sqrt(2),eps,f1,f2,recdepth+1,err);
		double Q2=Integrator24(f,(a+b)/2,b,acc/sqrt(2),eps,f3,f4,recdepth+1,err);
		return Q1+Q2;
	}
}

double Integrator( double f(double), double a, double b, double acc, double eps, double* err){
	double f2 = f(a+2*(b-a)/6);
	double f3 = f(a+3*(b-a)/6);
	int recdepth = 0;
	return Integrator24(f,a,b,acc,eps,f2,f3,recdepth,err);
}
