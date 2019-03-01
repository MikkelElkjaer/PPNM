#include<stdio.h>
#include<math.h>

double ln_over_sqrt(void);
double norm(double);
double H_integral(double); //make .h-file instead?

int main(){
	// problem 1
	fprintf(stderr,"Result of problem 1 is: %.18g \n \n \n", ln_over_sqrt());
	//printf("\nResult of problem 1 is: %.18g \n \n", ln_over_sqrt());

	//problem 2
	for(double alpha=0.05;alpha<3;alpha+=0.05){ //avoid alpha=0 as norm is then non-normalizable...
		double norm_value = norm(alpha);
		double H_int = H_integral(alpha);
		double E = H_int/norm_value;
		fprintf(stderr,"%g %g \n", alpha, E);
	}
return 0;
}
