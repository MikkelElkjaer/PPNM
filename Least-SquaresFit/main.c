#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_blas.h>

void lsfit(double (*f)(int, double), gsl_vector *, gsl_vector *, gsl_vector *, gsl_vector *, gsl_matrix *);

double fitfuncs(int i, double x){
	switch(i){
		case 0: return log(x); break;
		case 1: return 1.; break;
		case 2: return x; break;
		default: {fprintf(stderr,"%i is not a valid index for a fitting function. \n",i); return NAN;}
	}
}

int main(int argc, char** argv){
	int n = argc>1? atoi(argv[1]):0; //number of data points
	if(n == 0){printf("Error: main needs the number of data points as an argument \n"); return 0;}
	gsl_vector *x = gsl_vector_alloc(n), *y = gsl_vector_alloc(n), *dy = gsl_vector_alloc(n);
	for(int i=0;i<n;i++){
		double xi, yi, dyi; scanf("%lg %lg %lg", &xi, &yi, &dyi);
		gsl_vector_set(x,i, xi); gsl_vector_set(y,i,yi); gsl_vector_set(dy,i,dyi);
	}
	int m = 3; // number of fitting functions
	gsl_matrix * Covariance = gsl_matrix_alloc(m,m);
	gsl_vector * c = gsl_vector_alloc(m); // holds the coefficients of the fitting functions

	lsfit( &fitfuncs, x, y, dy, c, Covariance);
	printf("These are the fitting coefficients c0, c1, c2: \n");
	for(int i=0;i<m;i++) printf("%g ",gsl_vector_get(c,i));
	printf("\n \n");
	/* print also the fitting coefficients plus/minus their standard deviations: */
	printf("c0, c1, c2 plus/minus their standard deviation: \n");
	for(int i=0;i<m;i++) printf("%g ", gsl_vector_get(c,i) + sqrt(gsl_matrix_get(Covariance,i,i)));
	printf("\n");
	for(int i=0;i<m;i++) printf("%g ", gsl_vector_get(c,i) - sqrt(gsl_matrix_get(Covariance,i,i)));
        printf("\n \n");
	printf("The Covariance matrix: \n");
	for(int i=0;i<m;i++){
		for(int j=0;j<m;j++){
			printf("%8.5g ", gsl_matrix_get(Covariance,i,j));
		}
		printf("\n");
	}

	/* make output file for plotting the fit: */
	for(double z=gsl_vector_get(x,0);z<=gsl_vector_get(x,n-1);z+=0.2){
		double yfit = gsl_vector_get(c,0)*log(z) + gsl_vector_get(c,1) + gsl_vector_get(c,2)*z;
		fprintf(stderr,"%g %g ", z, yfit);
		for(int k=0;k<m;k++){
			double yfitp = 0, yfitm=0;
			for(int i=0;i<m;i++){
				double cp=gsl_vector_get(c,i);
				double cm=cp;
				if(i==k) cp+=sqrt(gsl_matrix_get(Covariance,i,i));
				if(i==k) cm-=sqrt(gsl_matrix_get(Covariance,i,i));
				yfitp+=fitfuncs(i,z)*cp;
				yfitm+=fitfuncs(i,z)*cm;
			}
			fprintf(stderr,"%g %g ",yfitp, yfitm);
		}
		fprintf(stderr,"\n");
	}
gsl_vector_free(x);gsl_vector_free(y);gsl_vector_free(dy);gsl_vector_free(c); gsl_matrix_free(Covariance);
return 0;
}
