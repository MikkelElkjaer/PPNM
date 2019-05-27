#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_blas.h>
#include<stdio.h>
#include<math.h>

int Newton(void f(gsl_vector*, gsl_vector*, gsl_matrix*), gsl_vector*, double);
int Quasi_Newton( void f(gsl_vector*, gsl_vector*), gsl_vector*, double);

void Rosenbrock(gsl_vector *, gsl_vector *, gsl_matrix *);
void Himmelblau(gsl_vector *, gsl_vector *, gsl_matrix *);
void RosenbrockNoH(gsl_vector *, gsl_vector *);
void HimmelblauNoH(gsl_vector *, gsl_vector *);

int main(){
        const int dim = 2; // dimension of the problems
        double eps = 1e-5;
        gsl_vector *x = gsl_vector_alloc(dim);

	/* Part A: Newton's method */
	printf("Newton's method with back-tracking linesearch:\n");
	gsl_vector_set(x,0,1.5); gsl_vector_set(x,1,1.5);
	int steps = Newton(Rosenbrock,x,eps);
	printf("Minimum of the Rosenbrock function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i, using precision: eps=%g, and starting point: (1.5,1.5) \n", steps, eps);

	/* minimum of Himmelblau function: (has one minimum at (3,2)) */
        gsl_vector_set(x,0,4.4); gsl_vector_set(x,1,4.2);
        steps = Newton(Himmelblau,x,eps);
        printf("Minimum of the Himmelblau function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i, using precision: eps=%g, and starting point: (4.4,4.2) \n", steps, eps);


	/* Part B: Quasi-Newton method with Broyden's update */
	printf("\n\nQuasi-Newton method with Broyden's update: \n(Note that the starting points have to be closer to the minima!)\n");
	gsl_vector_set(x,0,1.2); gsl_vector_set(x,1,1.2);
        steps = Quasi_Newton(RosenbrockNoH,x,eps);
        printf("Minimum of the Rosenbrock function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i, using precision: eps=%g, and starting point: (1.2,1.2) \n", steps, eps);

        /* minimum of Himmelblau function: (has one minimum at (3,2) - this particular guess yields another though: (-3.78,-3.28) */
        gsl_vector_set(x,0,4); gsl_vector_set(x,1,3.8);
        steps = Quasi_Newton(HimmelblauNoH,x,eps);
        printf("Minimum of the Himmelblau function is at: x=%.8g, y=%.8g\n",
                gsl_vector_get(x,0), gsl_vector_get(x,1));
        printf("Number of steps: %i, using precision: eps=%g, and starting point: (4,3.8) \n", steps, eps);

// HUSK FREE
gsl_vector_free(x);
return 0;
}
