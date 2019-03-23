#include<gsl/gsl_odeiv2.h>
#include<gsl/gsl_errno.h>
#include<math.h>

int ode_errfct(double x, const double u[], double dudx[], void* params){

        dudx[0] = 2/sqrt(M_PI) * exp(-x*x);

return GSL_SUCCESS;
}

double my_errfct(double x){
        gsl_odeiv2_system sys;
        sys.function = ode_errfct;
        sys.jacobian = NULL;
        sys.dimension = 1;
        sys.params = NULL;

        gsl_odeiv2_driver* driver;
        double hstart = 0.01, eps = 1e-6, abs = 1e-6;
	if(x<0) hstart*=-1;
        driver = gsl_odeiv2_driver_alloc_y_new(&sys,gsl_odeiv2_step_rkf45,hstart,abs,eps);

        //initial conditions:
        double x0=0;
        double u[] = {0};

        gsl_odeiv2_driver_apply(driver,&x0,x,u);

gsl_odeiv2_driver_free(driver);
return u[0];
}

