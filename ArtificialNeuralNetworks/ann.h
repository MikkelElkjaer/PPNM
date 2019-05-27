#include<gsl/gsl_vector.h>
#ifndef NEURONS_MADE
#define NEURONS_MADE // nice little trick to make this header file idempotent
typedef struct { int n; double (*f)(double); gsl_vector* data;} ann;
// data is kept in this form: data = {a_0,b_0,w_0,a_1,b_1,w_1,a_2,...}

// the following functions are found in ann.c
ann* ann_alloc(int n, double(*f)(double));
void ann_free(ann* network);
double ann_feed_forward(ann* network, double x);
void ann_train(ann* network, gsl_vector* xlist, gsl_vector* ylist);
#endif
