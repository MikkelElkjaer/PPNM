#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdbool.h>

struct counters {int steps; int fcalls;};
void QR_GS_decomp(gsl_matrix *, gsl_matrix *);
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);

struct counters Newton_with_Jacobian(void f(gsl_vector *x, gsl_vector *fx, gsl_matrix *J), gsl_vector *x, double eps){
/* Returns the approximate root of a nonlinear system of equations to a precision of eps using a user provided Jacobian. */
/* The inputs are: A function defining the system and calculating the Jacobian. A vector x with an initial guess which has
to be somewhat close to the actual root. This vector holds the approximate root upon completion. The error tolerance eps. */
	const int m = x->size;
        gsl_vector *fx = gsl_vector_alloc(m);
        gsl_matrix *J = gsl_matrix_alloc(m,m);
	gsl_vector *Dx = gsl_vector_alloc(m); // Delta x
	gsl_vector *xnew = gsl_vector_alloc(m);
	gsl_vector *fxnew = gsl_vector_alloc(m);
	double fxnorm, fxnewnorm;

	f(x,fx,J); // initialize function values and the Jacobian
	gsl_matrix *R = gsl_matrix_alloc(m,m); // upper triangular matrix for the QR-decomp routine
	int countsteps = 0;
	int functioncalls = 1;
	while(true){ countsteps++;
		QR_GS_decomp(J,R); // now J becomes Q in J=Q*R
		QR_GS_solve(J,R,fx,Dx);
		gsl_vector_scale(Dx,-1); // as we used fx in the solver above and not the actual -fx

		/* Do linesearch: */
		double lambda = 2;
		gsl_blas_ddot(fx,fx,&fxnorm);

		while(lambda>0.01){ functioncalls++;
			lambda/=2;
			for(int i=0;i<m;i++) gsl_vector_set(xnew,i, gsl_vector_get(x,i)+gsl_vector_get(Dx,i)*lambda);
			f(xnew,fxnew,J);
			gsl_blas_ddot(fxnew,fxnew,&fxnewnorm);
			if(fxnewnorm < (1-lambda/2)*fxnorm) break;
		}
		for(int i=0;i<m;i++){
			gsl_vector_set(x,i, gsl_vector_get(xnew,i));
			gsl_vector_set(fx,i, gsl_vector_get(fxnew,i));
		}
		if(fxnewnorm < eps || countsteps > 1000) break; // is 1000 a bit much maybe?
	}
struct counters c = {countsteps, functioncalls};
gsl_vector_free(fx);gsl_vector_free(Dx);gsl_vector_free(xnew);gsl_vector_free(fxnew);gsl_matrix_free(J);gsl_matrix_free(R);
return c;
}
