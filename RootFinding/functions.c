#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>

void system1(gsl_vector *x, gsl_vector *fx, gsl_matrix *J){
        double A = 10000;
        double y = gsl_vector_get(x,1), z = gsl_vector_get(x,0); // renaming x->z with respect to the exercise...
        gsl_vector_set(fx,0, A*z*y-1);
        gsl_vector_set(fx,1, exp(-z)+exp(-y)-1-1/A);
        gsl_matrix_set(J,0,0, A*y); gsl_matrix_set(J,0,1, A*z);
        gsl_matrix_set(J,1,0, -exp(-z)); gsl_matrix_set(J,1,1, -exp(-y));
return;}
void Rosenbrock(gsl_vector *x, gsl_vector *fx, gsl_matrix *J){
        double b = 100;
        double y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, -2*(1-z)-2*2*b*z*(y-z*z));
        gsl_vector_set(fx,1, 2*b*(y-z*z));
        gsl_matrix_set(J,0,0, 2-4*b*(y-z*z - 2*z*z));
        gsl_matrix_set(J,0,1, -4*b*z);
        gsl_matrix_set(J,1,0, -4*b*z);
        gsl_matrix_set(J,1,1, 2*b);
return;}
void Himmelblau(gsl_vector *x, gsl_vector *fx, gsl_matrix *J){
        double a = 11, b = 7, y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, 2*(z*z+y-a)*2*z + 2*(z+y*y-b));
        gsl_vector_set(fx,1, 2*(z*z+y-a) + 2*(z+y*y-b)*2*y);
        gsl_matrix_set(J,0,0, 4*(z*z+y-a) + 4*z*2*z + 2);
        gsl_matrix_set(J,0,1, 4*z + 4*y);
        gsl_matrix_set(J,1,0, 4*z + 4*y);
        gsl_matrix_set(J,1,1, 2 + 4*(z+y*y-b) + 4*y*2*y);
return;}

void system1NoJ(gsl_vector *x, gsl_vector *fx){
        double A = 10000;
        double y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, A*z*y-1);
        gsl_vector_set(fx,1, exp(-z)+exp(-y)-1-1/A);
return;}
void RosenbrockNoJ(gsl_vector *x, gsl_vector *fx){
        double b = 100;
        double y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, -2*(1-z)-2*2*b*z*(y-z*z));
        gsl_vector_set(fx,1, 2*b*(y-z*z));
return;}
void HimmelblauNoJ(gsl_vector *x, gsl_vector *fx){
        double a = 11, b = 7, y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, 2*(z*z+y-a)*2*z + 2*(z+y*y-b));
        gsl_vector_set(fx,1, 2*(z*z+y-a) + 2*(z+y*y-b)*2*y);
return;}

