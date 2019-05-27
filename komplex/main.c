#include"komplex.h"
#include"stdio.h"
#define TINY 1e-6

int main(){
	komplex a = {1,2}, b = {3,4};

	printf("testing komplex_add...\n");
	komplex r = komplex_add(a,b);
	komplex R = {4,6};
	komplex_print("a=",a);
	komplex_print("b=",b);
	komplex_print("a+b should   = ", R);
	komplex_print("a+b actually = ", r);

/* the following is optional */

	if( komplex_equal(R,r,TINY,TINY) )
		printf("test 'add' passed :) \n");
	else
		printf("test 'add' failed: debug me, please... \n");


// test own functions
	double acc = 10;
	double eps = 10;
	printf("testing komplex_equal... \n");
	int l = komplex_equal(a,b,acc,eps);
	printf("the return from komplex_equal was %i \n", l);

	printf("testing komplex_abs... \n");
	double size = komplex_abs(a);
	printf("abs of %g+i*%g is %g \n",a.re,a.im, size);

	printf("testing komplex_conjugate... \n");
	printf("The cc of %g+i*(%g) is %g+i*(%g) \n",a.re,a.im, komplex_conjugate(a).re, komplex_conjugate(a).im);
}
