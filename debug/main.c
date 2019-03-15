#include"stdio.h"
#include<gsl/gsl_matrix.h>

int print_half_00(gsl_matrix* m)
{
	double half = 1./2;
	int status = printf( "half m_{00} = %i\n", (int) (gsl_matrix_get(m,0,0)*half) );
	return status;
}

int main(void)
{
	gsl_matrix* m = gsl_matrix_alloc(1,1);
	gsl_matrix_set(m,0,0,66);
	printf("half m_{00} (should be 33):\n");
	int status = print_half_00(m);
	if(status==0)
		printf("status=%i : SOMETHING WENT TERRIBLY WRONG (status=0)\n",status);
	else
		printf("status=%i : everything went just fine (status>0)\n",status); //as far as I can tell, printf returning a 0 would correspond
// to nothing being printed, so I believe the if-else conditions had to be swapped...
	gsl_matrix_free(m);
return 0;
}
