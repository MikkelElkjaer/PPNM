#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>


int ode_schroe(double x, const double y[], double dydx[], void *params)
{
	double e = *(double*)params;
	dydx[0] = y[1];
	dydx[1] = -2*(e+ 1/x)*y[0];
	return GSL_SUCCESS;
}

double schroe_sol(double e,double xmax)
{
	gsl_odeiv2_system sys;
	sys.function = ode_schroe;
	sys.jacobian = NULL;
	sys.dimension = 2;
	sys.params = (void*)&e;

	gsl_odeiv2_driver *driver;
	double hstart = 0.1, abs = 1e-6, eps = 1e-6;
	driver = gsl_odeiv2_driver_alloc_y_new(&sys,
					       gsl_odeiv2_step_rkf45,
					       hstart, abs, eps);

	double rmin = 0.001; /* plateau of stabile solution seems to be around rmin = 0.001 with at least an order of magnitude in rmin*/
	double y[] = { rmin-rmin*rmin, 1-2*rmin }; /* boundary conditions for the wave function and its derivative at rmin */
	gsl_odeiv2_driver_apply(driver, &rmin, xmax, y);

	gsl_odeiv2_driver_free(driver);
	return y[0];
}
