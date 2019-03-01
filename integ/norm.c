#include<math.h>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_errno.h>

double integrand_norm(double x, void* params){
        double alpha = *(double *) params;
        double psi_normsquare = exp(-alpha*x*x);
        return psi_normsquare;
}

double norm(double alpha){
        int limit = 100;

        gsl_integration_workspace * w;
        w  = gsl_integration_workspace_alloc(limit);

        gsl_function F;
        F.function = &integrand_norm;
        F.params = (void*)&alpha;

        double result, error, acc=1e-7, eps=1e-7;
        int flag = gsl_integration_qagi(&F, acc, eps, limit, w, &result, &error);

        gsl_integration_workspace_free(w);

        if(flag!=GSL_SUCCESS) return NAN;
        return result;
}

