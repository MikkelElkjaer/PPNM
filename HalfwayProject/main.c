#include<stdio.h>
#include<math.h>

double myexp(double);

int main(void){
	for(double x=-10; x<10;x+=0.1){ /* quite a wide range but the result is still reasonably close to the exact exponentail fct. */
	printf("%g %g %g \n", x, myexp(x), exp(x));
	}
return 0;
}
