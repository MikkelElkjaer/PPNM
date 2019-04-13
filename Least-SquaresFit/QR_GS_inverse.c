#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);

void QR_GS_inverse(const gsl_matrix * Q, const gsl_matrix * R, gsl_matrix * B){
	// A (and thus Q) is assumed to be a square matrix in this problem
	gsl_vector * x = gsl_vector_alloc(Q->size2);
	gsl_vector * e_i = gsl_vector_alloc(Q->size2);
	gsl_matrix_set_identity(B);
	for(int i=0;i< Q->size2;i++){
		gsl_matrix_get_col(e_i,B,i);
		QR_GS_solve(Q, R, e_i, x);
		gsl_matrix_set_col(B,i,x);
	}
return;
}
