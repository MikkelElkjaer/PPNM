#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double integrand(double x, void* params){

	double f = log(x)/sqrt(x);
	return f;
}

double ln_over_sqrt(void){
	int limit = 100;

	gsl_integration_workspace * w;
        w  = gsl_integration_workspace_alloc(limit);

	gsl_function F;
        F.function = &integrand; // pointer in principle unneccessary
	F.params = NULL;

	double result, error, acc=1e-8, eps=1e-8;
        int flag = gsl_integration_qags (&F, 0, 1, acc, eps, limit, w, &result, &error);

        gsl_integration_workspace_free(w);

        if(flag!=GSL_SUCCESS) return NAN;
        return result;

}
