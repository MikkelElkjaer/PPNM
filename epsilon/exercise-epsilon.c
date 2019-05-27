#include<stdio.h>
#include<limits.h>
#include<float.h>

void part1(void);
void eps_part2(void);
int eps_part3();


int main()  {

part1();
eps_part2();

double a = 1.0;
double b = 1.01;
double tau = 0.05;
double epsilon = 0.005;

eps_part3(a, b, tau, epsilon);
return 0;
}

void part1(void) {
// exercise part i)
printf("INT_MAX is %i\n", INT_MAX);

int i=1;
while(i+1>i) i++;
printf("My max int is %i\n",i);

int counter = 1;
for(int j=1;j<INT_MAX; j++) counter++;
printf("INT_MAX in the for-loop reached: %i\n",counter);

int k=1;
do k++;  while(k<k+1) ;
printf("and then k reached: %i \n",k);

// part ii)
i=1;
printf("INT_MIN is %i \n", INT_MIN);
while(i-1<i) i--;
printf("My min int is %i \n",i);

counter = 1;
for(int j=1; j>INT_MIN;j--) counter--;
printf("INT_MIN in the for-loop reached: %i\n", counter);

k = 1;
do k--; while(k-1<k);
printf("and k reached: %i \n",k);

// part iii) - some iterative cases omitted...
float x = 1;
while(1+x!=1) {x/=2;} x*=2;
printf("FLT_EPS= %.10f, while loop gives: %.10f \n", FLT_EPSILON, x);
x=1; //and now using for-loop
for(x=1; 1+x!=1; x/=2) {} x*=2;
printf("for loop gives: %.10f \n",x);

double e = 1;
for(e=1; 1+e!=1; e/=2) {} e*=2;
printf("DBL_EPSILON = %.20g, while for-loop gives: %.20g \n", DBL_EPSILON, e);

long double t=1;
while(1+t!=1) {t/=2;} t*=2;
printf("LDBL_EPSILON = %.25Lg, while while-loop gives: %.25Lg \n", LDBL_EPSILON, t);
}

