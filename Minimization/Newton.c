#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdbool.h>

void QR_GS_decomp(gsl_matrix *, gsl_matrix *);
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);

int Newton(void f(gsl_vector *x, gsl_vector *df, gsl_matrix *H), gsl_vector *x, double eps){
/* Returns the coordinates of a minimum of the provided function. The user must provide the equations for the Hessian matrix.
The inputs are: A function defining the system and calculating the derivatives as well as the Hessian matrix.
A vector x with an initial guess which has to be somewhat close to the actual root.
This vector holds the approximate root upon completion. The error tolerance eps. */
	const int m = x->size;
        gsl_vector *df = gsl_vector_alloc(m);
        gsl_matrix *H = gsl_matrix_alloc(m,m);
	gsl_vector *Dx = gsl_vector_alloc(m); // Delta x
	gsl_vector *xnew = gsl_vector_alloc(m);
	gsl_vector *dfnew = gsl_vector_alloc(m);
	double dfnorm, dfnewnorm;
	double alpha = 0.33; // how strict of a Armijo condition to be used - could in principle be a parameter to the routine...
	double Armijoterm;

	f(x,df,H); // initialize the derivatives and the Hessian
	gsl_matrix *R = gsl_matrix_alloc(m,m); // upper triangular matrix for the QR-decomp routine
	int countsteps = 0;
	while(true){ countsteps++;
		QR_GS_decomp(H,R); // now H becomes Q in H=Q*R
		QR_GS_solve(H,R,df,Dx);
		gsl_vector_scale(Dx,-1); // as we used df in the solver above and not the actual -df

		/* Do linesearch: */
		double lambda = 2;
		gsl_blas_ddot(df,df,&dfnorm);

		while(lambda>0.01){
			lambda/=2;
			gsl_vector_scale(Dx,lambda);
			for(int i=0;i<m;i++) gsl_vector_set(xnew,i, gsl_vector_get(x,i)+gsl_vector_get(Dx,i));
			f(xnew,dfnew,H); // maybe implement a switch into the function to determine whether H should be updated or not...
			gsl_blas_ddot(dfnew,dfnew,&dfnewnorm);
			gsl_blas_ddot(Dx,df,&Armijoterm);
			gsl_vector_scale(Dx,1/lambda);
			if(dfnewnorm < dfnorm + alpha*Armijoterm) break;
		}
		for(int i=0;i<m;i++){
			gsl_vector_set(x,i, gsl_vector_get(xnew,i));
			gsl_vector_set(df,i, gsl_vector_get(dfnew,i));
		}
		if(dfnewnorm < eps || countsteps > 1000) break;
	}
gsl_vector_free(df);gsl_vector_free(Dx);gsl_vector_free(xnew);gsl_vector_free(dfnew);gsl_matrix_free(H);gsl_matrix_free(R);
return countsteps;
}
