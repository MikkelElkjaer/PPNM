#include<stdio.h>
#include"komplex.h"
#include<math.h> //nødvendig nu hvor tgmath er included?
#include<tgmath.h>

void komplex_print (char *s, komplex a) {
	printf ("%s (%g,%g)\n", s, a.re, a.im);
}

komplex komplex_new (double x, double y) {
	komplex z = { x, y };
	return z;
}

void komplex_set (komplex* z, double x, double y) {
	(*z).re = x;
	(*z).im = y;
}

komplex komplex_add (komplex a, komplex b) {
	komplex result = { a.re + b.re , a.im + b.im };
	return result;}

komplex komplex_sub (komplex a, komplex b) {
	komplex result = {a.re - b.re , a.im - b.im};
	return result;}

int komplex_equal (komplex a, komplex b, double acc, double eps) {
	double dist = sqrt( pow((a.re - b.re),2) + pow((a.im - b.im),2));
	if(dist<acc) {
	printf("the absolute difference is less than %g \n", acc); return 1;}
	else if( dist/( sqrt(pow(a.re,2)+pow(a.im,2)) + sqrt(pow(b.re,2)+pow(b.im,2))) < eps) {
	printf("the relative difference is less than %g \n", eps/2); return 1;}
	else return 0;
}

komplex komplex_conjugate (komplex z) {
	komplex result = {z.re , (-1)*z.im};
	return result;}

double  komplex_abs (komplex z) {
	double result = sqrt( pow(z.re,2) + pow(z.im,2));
	return result;}


