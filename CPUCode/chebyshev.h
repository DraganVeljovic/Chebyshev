/*
 * chebyshev.h
 *
 *  Created on: Sep 26, 2014
 *      Author: demo
 */

#ifndef CHEBYSHEV_H_
#define CHEBYSHEV_H_

struct Chebyshev {

	int n, m;

	double* c;

	double a, b;

	double (*fptr)(double);

};

void create (struct Chebyshev* c, double (*func)(double x), double aa, double bb, int nn);

void init (struct Chebyshev* c);

double eval (struct Chebyshev* c, double x, int m);

double* eval_vector(struct Chebyshev* c, double* x, int xSize, int m);
double* eval_vector_maxeler (struct Chebyshev* c, double* x, int xSize, int m);

void delete(struct Chebyshev* c);


#endif /* CHEBYSHEV_H_ */
