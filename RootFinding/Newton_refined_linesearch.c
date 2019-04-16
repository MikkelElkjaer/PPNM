#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<stdbool.h>

struct counters {int steps; int fcalls;};
void QR_GS_decomp(gsl_matrix *, gsl_matrix *);
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);

struct counters Newton_refined_linesearch( void f(gsl_vector* x, gsl_vector* fx), gsl_vector* x, double dx, double eps){
/* Returns the approximate root of a nonlinear system of equations to a precision of eps using a numerical estimate of the Jacobian.
The inputs are: A function defining the system. A vector x with an initial guess which has to be somewhat close to the actual root.
This vector holds the approximate root upon completion.
The finite difference value dx to be used in numerical evaluation of the Jacobian. The error tolerance eps. */
        const int m = x->size;
        gsl_vector *fx = gsl_vector_alloc(m);
        gsl_matrix *J = gsl_matrix_alloc(m,m);
        gsl_vector *Dx = gsl_vector_alloc(m); // Delta x
        gsl_vector *xnew = gsl_vector_alloc(m);
        gsl_vector *fxnew = gsl_vector_alloc(m);
        double fxnorm, fxnewnorm;

        f(x,fx); // initialize function values
        gsl_matrix *R = gsl_matrix_alloc(m,m); // upper triangular matrix for the QR-decomp routine
        int countsteps = 0;
        int functioncalls = 1;
        while(true){ countsteps++;
		for(int i=0;i<m;i++){ functioncalls++; // numerically calculate the Jacobian
			gsl_vector_set(x,i, gsl_vector_get(x,i)+dx);
			f(x,fxnew); //fxnew is temporarily used as a suiting container in calculating the Jacobian
			gsl_vector_sub(fxnew,fx); // fxnew is now a container for df (numerator in derivative estimate)
			for(int j=0;j<m;j++) gsl_matrix_set(J,j,i, gsl_vector_get(fxnew,j)/dx);
			gsl_vector_set(x,i, gsl_vector_get(x,i)-dx);
		}

                QR_GS_decomp(J,R); // now J becomes Q in J=Q*R
                QR_GS_solve(J,R,fx,Dx);
                gsl_vector_scale(Dx,-1); // as we used fx in the solver above and not the actual -fx

                /* Do refined linesearch using quadratic interpolation: */
                double lambda = 1;
                gsl_blas_ddot(fx,fx,&fxnorm);
		double g0 = 0.5*fxnorm;
		double g0prime = -fxnorm;

                while(lambda>0.005){ functioncalls++;
                        for(int i=0;i<m;i++) gsl_vector_set(xnew,i, gsl_vector_get(x,i)+gsl_vector_get(Dx,i)*lambda);
                        f(xnew,fxnew);
                        gsl_blas_ddot(fxnew,fxnew,&fxnewnorm);
			double gtrial = 0.5*fxnewnorm;
			double lambdatrial = lambda;
			double c = 1./(lambdatrial*lambdatrial)*(gtrial - g0 - g0prime*lambdatrial);
			lambda = -g0prime/(2*c);
                        if(fxnewnorm < (1-lambdatrial/2)*fxnorm) break;
                }
                for(int i=0;i<m;i++){ // update x and fx
                        gsl_vector_set(x,i, gsl_vector_get(xnew,i));
                        gsl_vector_set(fx,i, gsl_vector_get(fxnew,i));
                }
		double Dxnorm; gsl_blas_ddot(Dx,Dx,&Dxnorm);
                if(fxnewnorm < eps || Dxnorm < dx) break;
		if(countsteps > 1000){printf("\nError: Too many steps were taken. Please make a better initial guess. \n"); break;}
        }
struct counters c = {countsteps, functioncalls};
gsl_vector_free(fx);gsl_vector_free(Dx);gsl_vector_free(xnew);gsl_vector_free(fxnew);gsl_matrix_free(J);gsl_matrix_free(R);
return c;
}




