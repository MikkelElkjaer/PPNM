#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>

int ode_exp(double x, const double y[], double dydx[], void* params){

        dydx[0] = y[0];

return GSL_SUCCESS;
}

double myexp(double x){
       if (x == 0)
		return 1;
	if (x < 0)
		return 1 / myexp(-x);
	if (x > 2) {
		double z = myexp(x / 2);
		return z * z;
	}

	gsl_odeiv2_system sys;
        sys.function = ode_exp;
        sys.jacobian = NULL;
        sys.dimension = 1;
        sys.params = NULL;

        gsl_odeiv2_driver* driver;
        double hstart = 0.1, eps = 1e-5, abs = 1e-5;
        driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abs,eps);

        //initial conditions:
        double x0=0;
        double y[] = {1.0};

        gsl_odeiv2_driver_apply(driver,&x0,x,y);

gsl_odeiv2_driver_free(driver);
return y[0];
}

