#include<gsl/gsl_blas.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>

void QR_GS_decomp(gsl_matrix *, gsl_matrix *);
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);
void QR_GS_inverse(const gsl_matrix *, const gsl_matrix *, gsl_matrix *);

void lsfit( double (*fitfuncs)(int, double), gsl_vector * x, gsl_vector * y, gsl_vector * dy, gsl_vector * c, gsl_matrix * Covariance){
	int n = x->size; // number of data points
	int m = c->size; // number of fitting parameters
        gsl_matrix * A = gsl_matrix_alloc(n,m);
        gsl_vector * b = gsl_vector_alloc(n);
        for(int i=0;i<n;i++){ // initialize A and b to the fitting problem
                gsl_vector_set(b,i,gsl_vector_get(y,i)/gsl_vector_get(dy,i));
                for(int j=0;j<m;j++) gsl_matrix_set(A,i,j, (*fitfuncs)(j,gsl_vector_get(x,i))/gsl_vector_get(dy,i));
        }

        gsl_matrix * R = gsl_matrix_alloc(m,m);
	gsl_matrix * R_inv = gsl_matrix_alloc(m,m);
	gsl_matrix * Id = gsl_matrix_alloc(m,m); gsl_matrix_set_identity(Id);

        QR_GS_decomp(A,R); // turns A into Q in Q*R=A
        QR_GS_solve(A,R,b,c); // solves R*c = Q^T*b by direct back-substitution into c

	QR_GS_inverse(Id, R, R_inv); // calculates the inverse of R
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,R_inv,R_inv,0,Covariance); // calculates the covariance matrix

gsl_vector_free(b); gsl_matrix_free(A); gsl_matrix_free(R); gsl_matrix_free(Id); gsl_matrix_free(R_inv);
return;
}
