#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<pthread.h>

int Jacobi_sweep(gsl_matrix *, gsl_vector *, gsl_matrix *);

void print_matrix(gsl_matrix * M){
        const int n = M->size1, m = M->size2;
        for(int i=0;i<n;i++){
                for(int j=0;j<m;j++) printf("%4.3g ",gsl_matrix_get(M,i,j));
                printf("\n");
        }
return;
}

int main(int argc, char** argv){
	/* generate a symmetric, real matrix A: */
	int n = argc>1? atof(argv[1]):5;
	unsigned int seed = time(NULL);
	gsl_matrix *A = gsl_matrix_alloc(n,n);
	for(int i=0;i<n;i++){
		for(int j=i;j<n;j++){
		double a = (double)rand_r(&seed)/RAND_MAX*100;
		gsl_matrix_set(A,i,j,a);
		gsl_matrix_set(A,j,i,a);
		}
	}

	gsl_vector *eigvals = gsl_vector_alloc(n);
	gsl_matrix *V = gsl_matrix_alloc(n,n);
	int sweeps = Jacobi_sweep(A, eigvals, V);
	printf("The number of sweeps needed to reach converge for a particular %ix%i matrix is: %i \n",n,n,sweeps);

/* rebuild A from its lower triangular part: */
	for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) gsl_matrix_set(A,i,j, gsl_matrix_get(A,j,i));

/* check that V^T*A*V=D, where D is the diagonal matrix holding the eigenvalues in the eigvals vector: */
	gsl_matrix * D_check = gsl_matrix_alloc(n,n);
	gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,A,0,D_check);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,D_check,V,0,A); // save what is hopefully D in A...
	printf("This is V^T*A*V: \n");
	print_matrix(A);
	printf("The diagonal elements should match these: \n");
	for(int i=0;i<n;i++) printf("%.3g \n", gsl_vector_get(eigvals,i));

gsl_vector_free(eigvals); gsl_matrix_free(A); gsl_matrix_free(V); gsl_matrix_free(D_check);
return 0;
}
