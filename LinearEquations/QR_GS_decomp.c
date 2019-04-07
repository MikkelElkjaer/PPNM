#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>
#include<gsl/gsl_blas.h>
void QR_GS_decomp(gsl_matrix * A, gsl_matrix * R){
	int n = A->size1, m = A->size2; /* A is an (n x m) matrix, R is an (m x m) matrix */
	double R_ii, R_ij;
	gsl_vector * a_i = gsl_vector_alloc(n);
	gsl_vector * a_j = gsl_vector_alloc(n);
	for(int i=0;i<m;i++){
		gsl_matrix_get_col(a_i, A, i); // sets a_i to be the i'th coloumn of A
		gsl_blas_ddot(a_i, a_i, &R_ii); // calculates the dot product
		gsl_matrix_set(R,i,i, sqrt(R_ii)); // sets the i'th diagonal of R
		gsl_vector_scale(a_i,1./sqrt(R_ii)); // a_i -> q_i by scaling
		gsl_matrix_set_col(A,i, a_i); // stores the i'th coloumn of Q as the i'th coloumn in A.

		for(int j=i+1;j<m;j++){
			gsl_matrix_get_col(a_j,A,j);
			gsl_blas_ddot(a_i,a_j,&R_ij); // note that here a_i is in fact q_i due to the above scaling
			gsl_matrix_set(R,i,j,R_ij);
			gsl_vector_scale(a_i,R_ij); // scale q_i -> q_i*R_ij
			gsl_vector_sub(a_j,a_i); // update a_j
			gsl_matrix_set_col(A,j,a_j); // send the update into the matrix A aswell
			gsl_vector_scale(a_i,1./R_ij); // rescale q_i*R_ij -> q_i, so it may be used again for the next j
		}
	}
gsl_vector_free(a_i); gsl_vector_free(a_j);
return;
}
