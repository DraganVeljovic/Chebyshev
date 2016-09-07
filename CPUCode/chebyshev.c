/*
 * chebyshev.c
 *
 *  Created on: Sep 26, 2014
 *      Author: demo
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "chebyshev.h"

#include "Maxfiles.h"
#include "MaxSLiCInterface.h"

void create(struct Chebyshev* c, double(*func)(double x), double aa, double bb,
		int nn) {

	c->fptr = func;

	c->n = nn;
	c->m = nn;

	c->a = aa;
	c->b = bb;

	c->c = (double*) malloc(nn * sizeof(double));

}

void delete(struct Chebyshev* c) {
	free(c->c);
}

void init(struct Chebyshev* c) {

	int k, j;
	double fac, bpa, bma, y, sum;
	double *f;

	f = (double*) malloc(c->n * sizeof(double));

	bma = 0.5 * (c->b - c->a);
	bpa = 0.5 * (c->b + c->a);

	for (k = 0; k < c->n; k++) {
		y = cos(M_PI * (k + 0.5) / (double) c->n);
		f[k] = (*c->fptr)(y * bma + bpa);
	}

	fac = 2.0 / (double) c->n;

	for (j = 0; j < c->n; j++) {
		sum = 0.0;
		for (k = 0; k < c->n; k++)
			sum += f[k] * cos(M_PI * j * (k + 0.5) / (double) c->n);
		c->c[j] = fac * sum;

	}

	free(f);

}

double eval(struct Chebyshev* c, double x, int m) {

	double d = 0.0, dd = 0.0, y, y2;
	double sv;
	int j;

	// VAZNO
	//if ((x - a) * (x - b) > 0.0)
	//throw("x not in range in Chebyshev::eval");

	y2 = 2.0 * (y = (2.0 * x - c->a - c->b) / (c->b - c->a));

	for (j = m - 1; j > 0; j--) {
		sv = d;
		d = y2 * d - dd + c->c[j];
		dd = sv;
	}

	return y * d - dd + 0.5 * c->c[0];
}

double* eval_vector(struct Chebyshev* c, double* x, int xSize, int m) {

	double* returnValues = malloc(xSize * sizeof(double));

	for (int i = 0; i < xSize; i++) {
		//if ((x[i] - c->a) * (x[i] - c->b) > 0.0) return 0;
		returnValues[i] = eval(c, x[i], m);
	}

	return returnValues;
}

double* eval_vector_maxeler(struct Chebyshev* c, double* x, int xSize, int m) {

	double* returnValues = malloc(xSize * sizeof(double));
	double* coefToMap = malloc(m * sizeof(double));
	int i;

	for (i = 0; i < m; i++) {
		coefToMap[i] = c->c[m - 1 - i];
	}

	Chebyshev(c->a, c->b, m, xSize, x, returnValues, coefToMap);

	free(coefToMap);

	return returnValues;
}
