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

int next_exp(int lambda, int bound) {
	int x = -log(drand48()) / lambda;
	while (1) {
		if (ceil(x) > bound) {
			x = -log(drand48()) / lambda;
		} else {
			break;
		}
	}
	return x;
}

void assign(int pid, int** ID_TABLE) {
	for (int i = 0; i < ROWS; ++i) {
		for (int j = 0; j < COLS; ++j) {
			if (*(*(ID_TABLE+i)+j) == 0) {
				*(*(ID_TABLE+i)+j) = pid;
			}
		}
	}
}

int main(int argc, char** argv) {
	int n = atoi(*(argv+1));
	int n_CPU = atoi(*(argv+2));
	int seed = atoi(*(argv+3));
	srand48(seed);
	double lambda = atof(*(argv+4));
	int bound = atoi(*(argv+5));

	// keep count of children
	int children = n;
	// 0 if CPU-bound, 1 if IO-bound
	bool type = 0;

	// two character code table
	// if entry == 0, no process for that two character code
	int** ID_TABLE = calloc(ROWS, sizeof(int*));
	for (int i = 0; i < ROWS; ++i) {
		*(ID_TABLE+i) = calloc(COLS, sizeof(int));
	}

	pid_t p = -1;
	for (int i = 0; i < n; ++i) {
		p = fork();

		if (p == -1) {
			perror("fork() failed");
			return EXIT_FAILURE;
		}

		if (i >= n_CPU-1) {
			type = 1;
		}

		if (p == 0) break;
	}

	if (p == 0) {	// CHILD
		//printf("CHILD: Assigning process id %d to ID table\n", getpid());
		assign(getpid(), ID_TABLE);

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
	} else {
		while (children > 0) {
			waitpid(-1, NULL, 0);
			children--;
		}
		printf("<<< PROJECT PART I\n");
		if (n_CPU == 1) {
			printf("<<< -- process set (n=%d) with 1 CPU-bound process\n", n);
		} else {
			printf("<<< -- process set (n=%d) with %d CPU-bound processes\n", n, n_CPU);
		}
		printf("<<< -- seed=%d; lambda=%.6f; bound=%d\n", seed, lambda, bound);
	}

	return EXIT_SUCCESS;
}