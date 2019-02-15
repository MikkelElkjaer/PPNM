#include<stdio.h>
#include"nvector.h"
#include<assert.h>

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

double nvector_dot_product(nvector* u, nvector* v){
	
}
