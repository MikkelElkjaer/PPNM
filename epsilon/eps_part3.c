#include<stdio.h>
#include<math.h>
int eps_part3(double a, double b, double tau, double epsilon)
{
if(fabs(a-b)<tau) {
	printf("abs(a-b)<tau \n"); return 1;}
else if((fabs(a-b)/(fabs(a)+fabs(b))) < epsilon/2) {
	double x = fabs(a-b)/(fabs(a)+fabs(b));
	printf("lhs is: %g \n rhs is: %g \n", x, epsilon/2);
	printf("relative precision case \n"); 
	return 1;}
else {
	printf("a and b are far apart \n"); return 0;}


}
