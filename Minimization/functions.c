#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>

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

void RosenbrockNoH(gsl_vector *x, gsl_vector *fx){
        double b = 100;
        double y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, -2*(1-z)-2*2*b*z*(y-z*z));
        gsl_vector_set(fx,1, 2*b*(y-z*z));
return;}
void HimmelblauNoH(gsl_vector *x, gsl_vector *fx){
        double a = 11, b = 7, y = gsl_vector_get(x,1), z = gsl_vector_get(x,0);
        gsl_vector_set(fx,0, 2*(z*z+y-a)*2*z + 2*(z+y*y-b));
        gsl_vector_set(fx,1, 2*(z*z+y-a) + 2*(z+y*y-b)*2*y);
return;}

