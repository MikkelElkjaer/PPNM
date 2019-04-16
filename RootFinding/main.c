#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>

struct counters {int steps; int fcalls;};
struct counters Newton_with_Jacobian( void f(gsl_vector*, gsl_vector*, gsl_matrix*), gsl_vector*, double);
struct counters Newton( void f(gsl_vector*, gsl_vector*), gsl_vector*, double, double);
struct counters Newton_refined_linesearch( void f(gsl_vector*, gsl_vector*), gsl_vector*, double, double);

void system1(gsl_vector *, gsl_vector *, gsl_matrix *);
void Rosenbrock(gsl_vector *, gsl_vector *, gsl_matrix *);
void Himmelblau(gsl_vector *, gsl_vector *, gsl_matrix *);
void system1NoJ(gsl_vector *, gsl_vector *);
void RosenbrockNoJ(gsl_vector *, gsl_vector *);
void HimmelblauNoJ(gsl_vector *, gsl_vector *);

int main(){
	const int dim = 2; // dimension of (incidentally all) problems
	double eps = 1e-5;
	gsl_vector *x = gsl_vector_alloc(dim);

	printf("Here is part A with analytical Jacobians:\n");
	/* system of equations 1: */
	gsl_vector_set(x,0,1); gsl_vector_set(x,1,0.005); // put start guess in x
	struct counters c1 = Newton_with_Jacobian(system1,x,eps);
	printf("Root of A*x*y=1, exp(-x)+exp(-y) = 1+1/A with A=10000 is: x=%.8g, y=%.8g\n",
		gsl_vector_get(x,0), gsl_vector_get(x,1));
	printf("Number of steps: %i \nNumber of function calls: %i \n", c1.steps, c1.fcalls);

	/* minimum of Rosenbrock function: (has minimum at (1,1)) */
	gsl_vector_set(x,0,0.5); gsl_vector_set(x,1,0.5);
	struct counters c2 = Newton_with_Jacobian(Rosenbrock,x,eps);
	printf("Minimum of the Rosenbrock function is at: x=%.8g, y=%.8g\n",
		gsl_vector_get(x,0), gsl_vector_get(x,1));
	printf("Number of steps: %i \nNumber of function calls: %i \n", c2.steps, c2.fcalls);

	/* minimum of Himmelblau function: (has minimum at (3,2)) */
	gsl_vector_set(x,0,4.4); gsl_vector_set(x,1,4.2);
	struct counters c3 = Newton_with_Jacobian(Himmelblau,x,eps);
	printf("Minimum of the Himmelblau function is at: x=%.8g, y=%.8g\n",
		gsl_vector_get(x,0), gsl_vector_get(x,1));
	printf("Number of steps: %i \nNumber of function calls: %i \n", c3.steps, c3.fcalls);


	/* Exercise part B with numerical Jacobians: */
	printf("\n\nHere is part B with numerical Jacobians:\n");
	double dx = 1e-4; // OBSERVATION: this is the largest dx-value where the routine exactly reproduces that of part A.
	gsl_vector_set(x,0,1); gsl_vector_set(x,1,0.005);
	struct counters c1NoJ = Newton(system1NoJ,x,dx,eps);
	printf("Root of A*x*y=1, exp(-x)+exp(-y) = 1+1/A with A=10000 is: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c1NoJ.steps, c1NoJ.fcalls);

	gsl_vector_set(x,0,0.5); gsl_vector_set(x,1,0.5);
        struct counters c2NoJ = Newton(RosenbrockNoJ,x,dx,eps);
        printf("Minimum of the Rosenbrock function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c2NoJ.steps, c2NoJ.fcalls);

	gsl_vector_set(x,0,4.4); gsl_vector_set(x,1,4.2);
        struct counters c3NoJ = Newton(HimmelblauNoJ,x,dx,eps);
        printf("Minimum of the Himmelblau function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c3NoJ.steps, c3NoJ.fcalls);

	/* Exercise part C with refined linesearch (quadratic interpolation) using numerical Jacobians: */
	printf("\n\nHere is part C with refined linesearch: \n");
	gsl_vector_set(x,0,1); gsl_vector_set(x,1,0.005);
        c1NoJ = Newton_refined_linesearch(system1NoJ,x,dx,eps);
        printf("Root of A*x*y=1, exp(-x)+exp(-y) = 1+1/A with A=10000 is: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c1NoJ.steps, c1NoJ.fcalls);

        gsl_vector_set(x,0,0.5); gsl_vector_set(x,1,0.5);
        c2NoJ = Newton_refined_linesearch(RosenbrockNoJ,x,dx,eps);
        printf("Minimum of the Rosenbrock function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c2NoJ.steps, c2NoJ.fcalls);

        gsl_vector_set(x,0,4.4); gsl_vector_set(x,1,4.2);
        c3NoJ = Newton_refined_linesearch(HimmelblauNoJ,x,dx,eps);
        printf("Minimum of the Himmelblau function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i \nNumber of function calls: %i \n", c3NoJ.steps, c3NoJ.fcalls);

gsl_vector_free(x);
return 0;
}
