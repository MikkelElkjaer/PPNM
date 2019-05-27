#include<stdio.h>
#include"nvector.h"
#include<assert.h>
#include<stdlib.h>

nvector* nvector_alloc(int n){
  nvector* v = malloc(sizeof(nvector));
  (*v).size = n;
  (*v).data = malloc(n*sizeof(double));
  if( v==NULL ) fprintf(stderr,"error in nvector_alloc\n");
  return v;
}

void nvector_free(nvector* v){
	free(v->data);
	free(v);
}

void nvector_set(nvector* v, int i, double value){
	assert(i>=0);
        assert(i< v->size);
	(*v).data[i]=value;
}

double nvector_get(nvector* v, int i){
	assert(i>=0);
        assert(i< v->size);
	return (*v).data[i];
}

int nvector_equal(nvector* u, nvector* v, double epsilon){
	assert(u->size == v->size);
	int result = 1; //1 if equal, 0 otherwise
	for(int i=0;i< v->size;i++) {
		if (double_equal(u->data[i] , v->data[i]), epsilon)
			continue;
		else{
			result = 0;
			break; }
	}
	return result;
}

double nvector_dot_product(nvector* u, nvector* v){
	assert(u->size == v->size);
	double result;
	for(int i=0;i<(*v).size;i++) {result += u->data[i] * v->data[i];}
	return result;
}

void nvector_add(nvector* u, nvector* v) {  //note that the first vector is updated rather than a new vector is being created!!
	assert(u->size == v->size);
	for(int i=0;i< u->size;i++) {u->data[i] += v->data[i];}
}
