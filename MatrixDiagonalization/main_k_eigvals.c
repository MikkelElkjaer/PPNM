#include<stdio.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdlib.h>
#include<pthread.h>
#include<time.h>

int Jacobi_k_eigvals(gsl_matrix *, gsl_vector *, gsl_matrix *, int, int);
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
/* May take up to three inputs: matrix-dimension (nxn), number of wanted eigenvalues and 1/-1 for ordering by smallest/largest eigenvalue */
        int n = argc>1? atof(argv[1]):5;
	int k = argc>2? atof(argv[2]):3;
	int direction = argc>3? atof(argv[3]):1;
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
	clock_t begin_k = clock();
        int sweeps = Jacobi_k_eigvals(A, eigvals, V, k, direction); // the last integer should be 1/-1 for lowest/highest eigenvalues
	clock_t end_k = clock();
	printf("This is from the k_eigenvalues-routine: \n");
	printf("This many sweeps where needed to get the %i smallest/largest eigenvalues: %i \n", k, sweeps);

/* rebuild A from its lower triangular part: */
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) gsl_matrix_set(A,i,j, gsl_matrix_get(A,j,i));

/* check that V^T*A*V=D, where D is the diagonal matrix holding the eigenvalues in the eigvals vector: */
        gsl_matrix * D_check = gsl_matrix_alloc(n,n);
	gsl_matrix * VTA = gsl_matrix_alloc(n,n);
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,A,0,VTA);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,VTA,V,0,D_check);
        printf("This is V^T*A*V: \n");
        print_matrix(D_check);
        printf("The first %i diagonal elements should match these: \n", k);
        for(int i=0;i<k;i++) printf("%.3g \n", gsl_vector_get(eigvals,i));

/* for comparison, the full-sweep method is performed for the same matrix A below: */
	clock_t begin_sweep = clock();
	sweeps = Jacobi_sweep(A, eigvals, V); // eigvals and V are reinitialized in the routine
	clock_t end_sweep = clock();
	printf("This is from the sweep-routine: \n");
	printf("This many sweeps where needed to get the full decomposition: %i \n", sweeps);
/* rebuild A from its lower triangular part: */
        for(int i=0;i<n;i++) for(int j=i+1;j<n;j++) gsl_matrix_set(A,i,j, gsl_matrix_get(A,j,i));
/* check that V^T*A*V=D, where D is the diagonal matrix holding the eigenvalues in the eigvals vector: */
        gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,V,A,0,VTA);
        gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,VTA,V,0,D_check);
        printf("This is V^T*A*V: \n");
        print_matrix(D_check);
        printf("The diagonal elements should match these: \n");
        for(int i=0;i<n;i++) printf("%.3g \n", gsl_vector_get(eigvals,i));

/* compare the time consumption of the two methods: */
fprintf(stderr,"value-by-value method: Square matrix of size %i. %i eigenvalues were found in %g seconds \n",
		n, k, (double)(end_k-begin_k)/CLOCKS_PER_SEC);
fprintf(stderr,"cyclic method: Square matrix of size %i. A full diagonalization was done in %g seconds \n",
		n, (double)(end_sweep-begin_sweep)/CLOCKS_PER_SEC);

gsl_vector_free(eigvals); gsl_matrix_free(A); gsl_matrix_free(V); gsl_matrix_free(VTA); gsl_matrix_free(D_check);

return 0;
}
