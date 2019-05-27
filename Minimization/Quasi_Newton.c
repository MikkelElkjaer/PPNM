#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdbool.h>

#include<stdio.h>
#include<math.h>

int Quasi_Newton(void f(gsl_vector *x, gsl_vector *df), gsl_vector *x, double eps){
/* The inputs are: A function defining the system and calculating the Jacobian. A vector x with an initial guess which has
to be somewhat close to the actual root. This vector holds the approximate root upon completion. The error tolerance eps. */
	const int m = x->size;
        gsl_vector *df = gsl_vector_alloc(m);
        gsl_matrix *B = gsl_matrix_alloc(m,m); // inverse Hessian
	gsl_vector *Dx = gsl_vector_alloc(m); // Delta x
	gsl_vector *xnew = gsl_vector_alloc(m);
	gsl_vector *dfnew = gsl_vector_alloc(m);
	gsl_vector *y = gsl_vector_alloc(m);
	gsl_vector *u = gsl_vector_alloc(m);
	double dfnorm, dfnewnorm;
	double alpha = 0.33; // how strict of a Armijo condition to be used - could in principle be a parameter to the routine...
	double Armijoterm;
	double Broyden; double divergecondition = 1e-6;

	f(x,df); // initialize the derivatives
	int countsteps = 0;
	//int functioncalls = 1;
	gsl_matrix_set_identity(B); // zeroth order approx
	while(true){ countsteps++;
		gsl_blas_dgemv(CblasNoTrans, -1, B, df, 0, Dx); // get Newton's step Dx

		/* Do linesearch: */
		double lambda = 2;
		gsl_blas_ddot(df,df,&dfnorm);
		while(lambda>0.001){ // smaller lambda allowed as the derivative of the Rosenbrock function is very large even at (1.5,1.5)
			lambda/=2;
			gsl_vector_scale(Dx,lambda);
			for(int i=0;i<m;i++) gsl_vector_set(xnew,i, gsl_vector_get(x,i)+gsl_vector_get(Dx,i));
			f(xnew,dfnew);
			gsl_blas_ddot(dfnew,dfnew,&dfnewnorm);
			gsl_blas_ddot(Dx,df,&Armijoterm);
			gsl_vector_scale(Dx,1./lambda);
			if(dfnewnorm < (dfnorm + alpha*Armijoterm)) break;
		}

		/* update the inverse Hessian matrix using Broyden's update: */
		for(int i=0;i<m;i++) gsl_vector_set(y,i, gsl_vector_get(dfnew,i) - gsl_vector_get(df,i));
		for(int i=0;i<m;i++) gsl_vector_set(u,i, gsl_vector_get(Dx,i)*lambda); // u is temporarily s wrt. the notes
		gsl_blas_ddot(u,y,&Broyden);
		if(abs(Broyden) > divergecondition){
			gsl_blas_dgemv(CblasNoTrans, -1, B, y, 1, u); // now u becomes the u from the notes
			gsl_vector_scale(u,1./Broyden); // u then plays the role of c from the notes
			for(int i=0;i<m;i++){ for(int j=0;j<m;j++){ // the actual update:
				gsl_matrix_set(B,i,j, gsl_matrix_get(B,i,j) + gsl_vector_get(u,i)*gsl_vector_get(Dx,j)*lambda);
			}}
		}
		else gsl_matrix_set_identity(B); // if the update starts to diverge we reset the inverse Hessian

		for(int i=0;i<m;i++){
			gsl_vector_set(x,i, gsl_vector_get(xnew,i));
			gsl_vector_set(df,i, gsl_vector_get(dfnew,i));
		}
		if(dfnewnorm < eps || countsteps > 10000) break;
	}
gsl_vector_free(df);gsl_vector_free(Dx);gsl_vector_free(xnew);
gsl_vector_free(dfnew);gsl_matrix_free(B);gsl_vector_free(u); gsl_vector_free(y);
return countsteps;
}
