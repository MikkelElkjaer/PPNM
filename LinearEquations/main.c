#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<stdlib.h>
#include<gsl/gsl_blas.h>
#include<pthread.h>
void QR_GS_decomp(gsl_matrix *, gsl_matrix *);
void QR_GS_solve(const gsl_matrix *, const gsl_matrix *, const gsl_vector *, gsl_vector *);
void QR_GS_inverse(const gsl_matrix *, const gsl_matrix *, gsl_matrix *);

gsl_matrix * get_random_matrix(int n, int m){
        unsigned int seed = time(NULL); /* make thread-local seed for rand_r()*/
	gsl_matrix * M = gsl_matrix_alloc(n,m);
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++) gsl_matrix_set(M,i,j,(double)rand_r(&seed)/RAND_MAX*100);
	}
return M;
}

void print_matrix(gsl_matrix * M){
	const int n = M->size1, m = M->size2;
	for(int i=0;i<n;i++){
		for(int j=0;j<m;j++) printf("%4.3g ",gsl_matrix_get(M,i,j));
		printf("\n");
	}
return;
}

int main(){
	/* part A - QR-decomposition: */
	int n = 5, m = 3;
	gsl_matrix * A = get_random_matrix(n,m);
	gsl_matrix * R = gsl_matrix_alloc(m,m);
	gsl_matrix_set_zero(R);

	gsl_matrix * A_copy = gsl_matrix_alloc(n,m);
	gsl_matrix_memcpy(A_copy,A);

	QR_GS_decomp(A,R); // turns A into Q and computes R; A = Q*R

	/* show that Q*R=A, Q^T*Q=1 and that R is upper diagonal: */
	gsl_matrix * QR = gsl_matrix_alloc(n,m);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A, R, 0, QR); // computes Q*R and puts it into QR (note: Q is called A)
	printf("This is Q*R: \n");
	print_matrix(QR);
	printf("This is A: \n");
	print_matrix(A_copy);

	gsl_matrix * QTQ = gsl_matrix_alloc(m,m);
        gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1, A, A, 0, QTQ); // use QTQ to store the matrix product of Q transposed and Q. Again, Q is called A
	printf("This is Q^T*Q: \n");
	print_matrix(QTQ);

	printf("This is R: \n");
	print_matrix(R);

	gsl_matrix_free(A); gsl_matrix_free(R); gsl_matrix_free(QR); gsl_matrix_free(A_copy); gsl_matrix_free(QTQ);

	/* solve Ax=b (A has to be square): */
	gsl_matrix * A2 = get_random_matrix(n,n);
        gsl_matrix * R2 = gsl_matrix_alloc(n,n);
        gsl_matrix_set_zero(R2);
	gsl_matrix * A2_copy = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(A2_copy,A2);

	QR_GS_decomp(A2,R2);
	gsl_vector *b = gsl_vector_alloc(n);
	gsl_vector *x = gsl_vector_alloc(n);
	unsigned int seed = time(NULL); /* make thread-local seed for rand_r()*/
	for(int i=0;i<n;i++) gsl_vector_set(b,i,(double)rand_r(&seed)/RAND_MAX*100);
	QR_GS_solve(A2,R2,b,x); // solves QRx=b. Still, Q is called A

	/* check that the solution x is correct: */
	gsl_vector *b_check = gsl_vector_alloc(n);
	gsl_blas_dgemv(CblasNoTrans, 1, A2_copy, x, 0, b_check); // compute A*x into b_check
	printf("\nSolving QR*x=b correctly should produce two identical vectors (b) here: \n");
	for(int i=0;i<n;i++) printf("%7.4g %7.4g \n", gsl_vector_get(b,i), gsl_vector_get(b_check,i));

	gsl_vector_free(b); gsl_vector_free(x); gsl_vector_free(b_check); gsl_matrix_free(A2); gsl_matrix_free(R2); gsl_matrix_free(A2_copy);
	// ALSO ASK: should you incorporate an assert to make the diag(R)>0?

	/* Part B (Matrix inverse by Gram-Schmidt QR-factorization): */
	gsl_matrix * A3 = get_random_matrix(n,n);
	gsl_matrix * R3 = gsl_matrix_alloc(n,n);
	gsl_matrix * B3 = gsl_matrix_alloc(n,n);
	gsl_matrix_set_zero(R3);
	gsl_matrix * A3_copy = gsl_matrix_alloc(n,n);
	gsl_matrix_memcpy(A3_copy,A3);

	QR_GS_decomp(A3,R3); // A3 is now Q in Q*R = A
	QR_GS_inverse(A3,R3,B3); // calculate the inverse and store it in B3

	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1, A3_copy, B3, 0, R3); // save A*A⁻¹ in R

	printf("\nThis is A*A⁻¹: \n");
	print_matrix(R3);

return 0;
}
