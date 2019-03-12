#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<pthread.h>

/* the 1x1 square is centered on (0,0). The inscribed circle has radius 0.5 and hitting within the circle is equivalent to
have a distance less than 0.5 to the origin. */
/* NB: "man pthreads" states that rand() is thread-safe... and "man rand_r()" states that rand_r() became obsolete in 2009.*/

struct count {int a, b ; int N_in;};

void* throw_darts(void* param){
	unsigned int seed = time(NULL); /* make thread-local seed for rand_r()*/
	struct count* data = (struct count*)param;
	int N_in = 0;
	int a = data->a;
	int b = data->b;
	double x, y;
	for(int i=a;i<b;i++){
		x = -0.5 + (double)rand_r(&seed) /(double)RAND_MAX;
		y = -0.5 + (double)rand_r(&seed) /(double)RAND_MAX;
		double d = sqrt(x*x+y*y);
		if (d <= 0.5) N_in++;
	}
	data->N_in = N_in;
}

int main(void){
	/* Using pthreads:*/
	for(int N=1e2; N<=1e8; N*=5){
	struct count lower, upper;
	lower.a = 0;
	lower.b = N/2;
	upper.a = N/2;
	upper.b = N;

	pthread_t t1;
	pthread_create(&t1,NULL,throw_darts,(void*)&lower);

	throw_darts((void*)&upper);

	pthread_join(t1,NULL);

	double my_pi = (double)4*(upper.N_in+lower.N_in)/N;
	printf("%3i %4.8lg %4lg\n",N,my_pi,M_PI-my_pi);
	}

	printf("\n \n");

	/* Using OpenMP: */
	for(int N=1e2; N<=1e8; N*=5){
        struct count lower, upper;
        lower.a = 0;
        lower.b = N/2;
        upper.a = N/2;
        upper.b = N;

		#pragma omp parallel sections
		{
		#pragma omp section
			{
        		throw_darts((void*)&lower);
			}
		#pragma omp section
			{
			throw_darts((void*)&upper);
			}
		}

        double my_pi = (double)4*(upper.N_in+lower.N_in)/N;
	printf("%3i %4.8lg %4lg\n",N,my_pi,M_PI-my_pi);
	}
return 0;
}
