#include<stdio.h>
#include<limits.h>

void eps_part2(void)
{
int max = INT_MAX/3;

float sum_up_float = 0.0;
int i;
for(i=1 ; i<max ; i++) {sum_up_float += 1.0f/i;}
printf("sum_up_float reached %.10f \n", sum_up_float);

float sum_down_float = 0.0;
for(i=0;i<max;i++) {sum_down_float += 1.0f/(max-i);}
printf("sum_down_float reached %.10f \n", sum_down_float);

// their difference is due to the fact that floats do not obey the associative rule.
// that is, when adding, say, 5e-12 to 1e10 and subtracting 1e10 afterwards the 5e-12 gets "washed out".

// part iii) does the sum converge as a fct of max? omitted...

// the same for doubles:
double sum_up_double = 0.0;
for(i=1 ; i<max ; i++) {sum_up_double += 1.0/i;}
printf("sum_up_double reached %.20g \n", sum_up_double);

double sum_down_double = 0.0;
for(i=0 ; i<max ; i++) {sum_down_double += 1.0/(max-i);}
printf("sum_down_double reached %.20g \n", sum_down_double);
printf("Note that the difference of the doubles is far less than for \n the floats (the washing out of precision is far less impactful \n");

}
