#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<math.h>

int Jacobi_k_eigvals(gsl_matrix *A, gsl_vector *eigvals, gsl_matrix *V, int k, int direction){
/* Jacobi eigenvalue decomposition of k eigenvalues. The direction should be 1/-1 for ordering by smallest/largest
eigenvalues. The vector eigvals stores the k eigenvalues and the matrix V stores the eigenvectors as
its coloumns. The upper triangular part of A is destroyed. The function returns the number of sweeps needed.
One may rebuild A from its lower triangular part utilizing that A is assumed symmetric and real. */
int n=A->size1;
int sweeps = 0;
int change_flag; // flags if a Jacobi-rotation results in a change of the eigenvalues affected
int i, p, q; // iterator and indicies

gsl_matrix_set_identity(V); // Maybe this should only be a block of dim k which is an eye:
//for(i=k;i<n;i++) gsl_matrix_set(V,i,i,0); // make zeros on all but the k first diagonal positions                  
for(i=0;i<n;i++) gsl_vector_set(eigvals,i, gsl_matrix_get(A,i,i)); // Easy way to keep track of the diagonal of A

// Also maybe truncate eigvals to the first k elements at the end??

for(p=0;p<k;p++){
do {
        change_flag = 0; sweeps++;
        //for(p=0;p<k;p++)
	for(q=p+1;q<n;q++) {
                double App = gsl_vector_get(eigvals,p);
                double Aqq = gsl_vector_get(eigvals,q);
                double Apq = gsl_matrix_get(A,p,q);
                double phi = 0.5*atan2(direction*2*Apq, direction*(Aqq-App)); // direction: 1/-1 for ordering by smallest/largest eigenvalues
                double c = cos(phi);
                double s = sin(phi);
                double App_new = c*c*App - 2*s*c*Apq + s*s*Aqq;
                double Aqq_new = s*s*App + 2*s*c*Apq + c*c*Aqq;
                if(App_new != App || Aqq_new != Aqq){ // update only if the Jacobi rotation changed the eigenvalues
                        change_flag = 1; // the flag is now true
                        gsl_vector_set(eigvals,p,App_new); // update the eigenvalues
                        gsl_vector_set(eigvals,q,Aqq_new);
                        gsl_matrix_set(A,p,q,0.0); // as a result of the rotation
                        /* update the rest of the upper part of A where one and only one index is p or q using eqn. (9) from the notes */
                        for(i=p+1;i<q;i++){
                                double Api = gsl_matrix_get(A,p,i); // elements horizontally between (p,p) and (p,q)
                                double Aiq = gsl_matrix_get(A,i,q); // elements vertically between (p,q) and (q,q)
                                gsl_matrix_set(A,p,i, c*Api - s*Aiq);
                                gsl_matrix_set(A,i,q, c*Aiq + s*Api);  }
                        for(i=q+1;i<n;i++){
                                double Api = gsl_matrix_get(A,p,i); // elements to the right of (p,q)
                                double Aqi = gsl_matrix_get(A,q,i); // elements to the right of (q,q)
                                gsl_matrix_set(A,p,i, c*Api - s*Aqi);
                                gsl_matrix_set(A,q,i, c*Aqi + s*Api);  }
                        /* update the eigenvectors using equation (11) from the notes */
                        for(i=0;i<n;i++){
                                double Vip = gsl_matrix_get(V,i,p);
                                double Viq = gsl_matrix_get(V,i,q);
                                gsl_matrix_set(V,i,p, c*Vip - s*Viq);
                                gsl_matrix_set(V,i,q, c*Viq + s*Vip);  }

                }
        }
} // closing "do"
while(change_flag != 0);
} // closing outer for
return sweeps;
}

