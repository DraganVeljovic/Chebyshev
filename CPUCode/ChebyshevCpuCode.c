#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "chebyshev.h"

double func(double x) {
	return cos(x);
}

int main(void)
{
	int max = 1;
	int min = -1;

	struct Chebyshev c;

	create(&c, &func, min, max, 256);

	init(&c);

	int randomDataAmount = 10000000;
	int precision = 256;
	
	double* randomData = malloc(randomDataAmount * sizeof(double));
	double *resultCPU, *resultMaxeler;
	
	for (int i = 0; i < randomDataAmount; i++) {
		randomData[i] = (max - min) * (rand() / (RAND_MAX * 1.0)) + min;
	}

	clock_t cpu_start = clock();

	resultCPU = eval_vector(&c, randomData, randomDataAmount, precision);

	clock_t cpu_end = clock();

	clock_t maxeler_start = clock();

	resultMaxeler = eval_vector_maxeler(&c, randomData, randomDataAmount, precision);

	clock_t maxeler_end = clock();

	printf("\nVreme CPU: %lf", (double)(cpu_end - cpu_start) / (double)CLOCKS_PER_SEC);
	printf("\nVreme Maxeler: %lf", (double)(maxeler_end - maxeler_start) / (double)CLOCKS_PER_SEC);

	for (int i = 0; i < randomDataAmount; i++) {
		if (abs(resultCPU[i] - resultMaxeler[i])) {
			printf("\n Neslaganje!");
			break;
		}
	}

	printf("\nGotovo!");

	delete(&c);

	free(randomData);
	free(resultCPU);
	free(resultMaxeler);

	return 0;
}
