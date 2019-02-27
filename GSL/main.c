#include<stdio.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_linalg.h>
#include<math.h>
void airy(void); // rather declare it here than make a .h-file for just one function
void matrix_print(const char*, gsl_matrix*);
void vector_print(const char*, gsl_vector*);

int main(){
	airy(); // make plot of airy functions

	// part II:
	int size1,size2;
	scanf("%i",&size1);
	scanf("%i",&size2);

	gsl_matrix* A = gsl_matrix_alloc(size1,size2);
	gsl_matrix_fscanf(stdin,A);

	matrix_print("A=",A);

	gsl_vector* b=gsl_vector_alloc(size2);
	gsl_vector* x = gsl_vector_alloc(size2);
	gsl_vector_fscanf(stdin,b);

	vector_print("b=",b);

	gsl_matrix* B = gsl_matrix_alloc(size1,size2); //as HH destroys A
	gsl_matrix_memcpy(B,A);
	gsl_linalg_HH_solve(A,b,x);
	vector_print("x=",x);

	gsl_blas_dgemv(CblasNoTrans,1,B,x,0,b);

	vector_print("Ax=",b); // if this line is true, the solution is correct

	gsl_matrix_free(A);
	gsl_matrix_free(B);
	gsl_vector_free(b);
	gsl_vector_free(x);

return 0;
}


void vector_print(const char* s, gsl_vector* v){
	printf("%s",s);
	for(int j=0;j< v->size;j++)
			printf("%8.3g ",gsl_vector_get(v,j));
	printf("\n");

}

void matrix_print(const char* s, gsl_matrix* A){
	printf("%s\n",s);
	for(int i=0; i<A->size1;i++){
		for(int j=0; j< A->size2;j++){
			printf("%8.3g ",gsl_matrix_get(A,i,j));
		}
		printf("\n");
	}
}
