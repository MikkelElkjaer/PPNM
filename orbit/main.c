#include<stdio.h>
#include<math.h>

double part1(double);
double part2(double, double*, double*);

int main(){
	//part1:
	for(double x=0;x<=3;x+=0.1){
	printf("%g %g \n", x, part1(x));
	}

	//part2:

	//part i)
	int n = 3; //number of initial conditions to be specified
	double init_conds[n];
	init_conds[0] = 0, init_conds[1] = 1, init_conds[2] = 0;
	double epsilon = 0;
	for(double x=0; x<2*M_PI+0.05;x+=0.05){
		fprintf(stderr,"%g %g \n", x, part2(x,init_conds,&epsilon));
	}
	fprintf(stderr,"\n \n");
	// part ii)
	init_conds[2] = -0.5;
	for(double x=0; x<2*M_PI+0.05;x+=0.05){
		fprintf(stderr,"%g %g \n", x, part2(x,init_conds,&epsilon));
	}
	fprintf(stderr,"\n \n");
	// part iii)
	epsilon = 0.025;
	for(double x=0; x<20*M_PI;x+=0.05){
		fprintf(stderr,"%g %g \n", x, part2(x,init_conds,&epsilon));
	}

return 0;
}
