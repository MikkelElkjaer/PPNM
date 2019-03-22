#include<gsl/gsl_multiroots.h>
#include<gsl/gsl_errno.h>
#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<gsl/gsl_vector.h>

double schroe_sol(double e, double rmax);

int root_equation
(const gsl_vector * x, void * params, gsl_vector * f)
{
	double rmax = *(double*)params;
	double e = gsl_vector_get(x,0);
	double mismatch = schroe_sol(e,rmax);
	gsl_vector_set(f,0,mismatch);
return GSL_SUCCESS;
}

int main(void){
	FILE* energystream = fopen("hydrogen_energies.txt","w");
	for(double rmax = 3; rmax<=8; rmax+=0.5) {
		double estart = -1;

		gsl_multiroot_function F;
		F.f=root_equation;
		F.n=1;
		F.params=(void*)&rmax;

		gsl_multiroot_fsolver * S;
		S = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids,F.n);

		gsl_vector* start = gsl_vector_alloc(F.n);
		gsl_vector_set(start,0,estart);
		gsl_multiroot_fsolver_set(S,&F,start);

		int flag;
		do{
			gsl_multiroot_fsolver_iterate(S);
			flag=gsl_multiroot_test_residual(S->f,1e-6);
		}while(flag==GSL_CONTINUE);

		double result=gsl_vector_get(S->x,0); /* that is, the resulting energy */
		gsl_multiroot_fsolver_free(S);
		gsl_vector_free(start);
		fprintf(energystream,"%g %g %g\n",rmax,result,schroe_sol(result,rmax)); /* rmax, energy, wavefunctions value at rmax */
		if(((int)(2*rmax))%4 == 0) { /* plot wavefct only for rmax = 4,6 and 8 */
			for (double r=0.001; r<=rmax;r+=rmax/100){
			fprintf(stderr,"%g %g %g,\n",r,schroe_sol(result,r),r*exp(-r)); /* to plot the wavefunction */
			};
                        fprintf(stderr,"\n \n"); /* make a new data set in hydrogen_data.txt */
		};
	};
	fclose(energystream);
return 0;
}
