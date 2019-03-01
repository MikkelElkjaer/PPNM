#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_part2(double t, const double y[], double dydt[], void* params){

	double epsilon = *(double *)params;

	dydt[0] = y[1];
	dydt[1] = 1-y[0]+epsilon*y[0]*y[0];

return GSL_SUCCESS;
}

double part2(double t, double* init_conds, double* epsilon){

	gsl_odeiv2_system sys;
	sys.function = ode_part2;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)epsilon;

	gsl_odeiv2_driver* driver;
	double hstart = 0.05, abs = 1e-6, eps = 1e-6;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abs,eps);

	//initial conditions:
	double t0 = init_conds[0];
	double y[2] = {init_conds[1] , init_conds[2]};

	gsl_odeiv2_driver_apply(driver,&t0,t,y);

	gsl_odeiv2_driver_free(driver);

return y[0];
}
