#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h> //necessary????

int ode_part1(double t, const double y[], double dydt[], void* params){

	dydt[0] = y[0]*(1-y[0]);

return GSL_SUCCESS; // in principle just 0 for success...
}


double part1(double t){
	gsl_odeiv2_system sys;
	sys.function = ode_part1;
	sys.jacobian = NULL;
	sys.dimension = 1;
	sys.params = NULL;

	gsl_odeiv2_driver* driver;
	double hstart = 0.1, eps = 1e-5, abs = 1e-5;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abs,eps);

	//initial conditions:
	double t0=0;
	double y[] = {0.5};

	gsl_odeiv2_driver_apply(driver,&t0,t,y);

gsl_odeiv2_driver_free(driver);
return y[0];
}
