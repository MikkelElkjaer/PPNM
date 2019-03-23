#include<stdio.h>
#include<stdlib.h>

double my_errfct(double);

int main(int argc, char** argv){
	double a = argc>1? atof(argv[1]):-2;
	double b = argc>2? atof(argv[2]):2;
	double dx = argc>3? atof(argv[3]):0.1;

	for(double x = a; x<=b; x+=dx)
		printf("%5.5g %5.5g \n", x, my_errfct(x));

return 0;
}
