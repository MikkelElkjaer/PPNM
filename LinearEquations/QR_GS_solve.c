#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>

void QR_GS_solve(const gsl_matrix *Q, const gsl_matrix *R, const gsl_vector *b, gsl_vector *x){
	gsl_blas_dgemv(CblasTrans, 1, Q, b, 0, x); // computes Q^T*b and saves it in x
	/* do in-place back-substitution: */
	for(int i=x->size -1;i>=0;i--){
		double a = gsl_vector_get(x,i);
		for(int j=i+1;j< x->size;j++) a -= gsl_matrix_get(R,i,j)*gsl_vector_get(x,j);
		gsl_vector_set(x,i,a/gsl_matrix_get(R,i,i));
	}
return;
}
