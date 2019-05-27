#include<stdio.h>
#include<math.h>
int double_equal(double a, double b, double epsilon)
{
if(fabs(a-b)<epsilon) {
	 return 1;}
else {
return 0;}
}
