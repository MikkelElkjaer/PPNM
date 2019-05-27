#include<stdio.h>
#include<math.h>
#include<complex.h>

int main() {
	double x = tgamma(5);
	double y = j1(0.5);
	complex a = csqrt(-2);
	complex b = cexp(I);
	complex c = cexp(I*M_PI);
	complex d = cpow(I,M_E);
	printf("x=%g, y=%g\n", x, y);
	printf("a=%g + I* %g\n", creal(a),cimag(a));
	printf("b=%g + I* %g \n", creal(b),cimag(b));
	printf("c=%g + I* %g \n", creal(c), cimag(c));
	printf("d=%g + I*%g \n", creal(d), cimag(d));

return 0;


}
