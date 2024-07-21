#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

#define ROWS 26
#define COLS 10

int next_exp(double lambda, double bound) {
	double r = drand48();
	double x = -log(r) / lambda;
	while (1) {
		if (ceil(x) > bound) {
			r = drand48();
			x = -log(r) / lambda;
		} else {
			break;
		}
	}
	return x;
}

// char* fetchCode(int pid, int** ID_TABLE) {
// 	int m = -1;
// 	int n = -1;
// 	for (int i = 0; i < ROWS; ++i) {
// 		for (int j = 0; j < COLS; ++j) {
// 			if (*(*(ID_TABLE+i)+j) == pid) {
// 				m = i;
// 				n = j;
// 				break;
// 			}
// 		}
// 	}
// 	char* code = calloc(2, sizeof(char));
// 	*code = m + '0'
// }

int main(int argc, char** argv) {
	if (argc != 6) {
		fprintf(stderr, "ERROR: invalid number of arguments\n");
		abort();
	}

	int n = atoi(*(argv+1));
	int n_CPU = atoi(*(argv+2));
	int seed = atof(*(argv+3));
	srand48(seed);
	double lambda = atof(*(argv+4));
	int bound = atof(*(argv+5));

	if (n <= 0) {
		fprintf(stderr, "ERROR: invalid number of processes\n");
		abort();
	}

	if (n_CPU < 0 || n_CPU > n) {
		fprintf(stderr, "ERROR: invalid number of processes that are CPU-bound\n");
		abort();
	}

	// if (seed) {
	// 	fprintf(stderr, "ERROR: invalid seed value\n");
	// 	abort();
	// }

	if (lambda <= 0) {
		fprintf(stderr, "ERROR: invalid lambda value\n");
		abort();
	}

	// if (bound) {
	// 	fprintf(stderr, "ERROR: invalid seed value\n");
	// 	abort();
	// }

	// keep count of children
	int children = n;
	// 0 if CPU-bound, 1 if IO-bound
	bool type = 0;

	printf("<<< PROJECT PART I\n");
	if (n_CPU == 1) {
		printf("<<< -- process set (n=%d) with 1 CPU-bound process\n", n);
	} else {
		printf("<<< -- process set (n=%d) with %d CPU-bound processes\n", n, n_CPU);
	}
	printf("<<< -- seed=%d; lambda=%.6f; bound=%d\n", seed, lambda, bound);

	for (int i = 0; i < n; ++i) {
		if (i >= n_CPU-1) {
			type = 1;
		}

		// Step 1
		int arrival_time = floor(next_exp(lambda, bound));

		// Step 2
		double uniform_distribution = drand48();
		int CPU_bursts = ceil(uniform_distribution) * 32;

		// Step 3
		for (int i = 0; i < CPU_bursts; ++i) {

			int CPU_burst_time = ceil(next_exp(lambda, bound));
			int IO_burst_time = -1;
			if (i == CPU_bursts-1) {	// do not generate IO burst time for the last CPU burst
				IO_burst_time = ceil(next_exp(lambda, bound));
			}
			if (type == 0) {	// CPU-bound
				CPU_burst_time *= 4;
			} else {			// IO-bound
				IO_burst_time *= 8;
			}
		}
	}

	return EXIT_SUCCESS;
}