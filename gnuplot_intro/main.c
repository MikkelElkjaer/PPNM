#include<stdio.h>
#include<math.h>
int main(void) {
	int n=20;
	double a =0.2, b=2*M_PI;
	double dx=(b-a)/(n-1);
	for( int i=0;i<n;i++) printf("%g %g \n",i*dx,sin(i*dx));
return 0;
}
