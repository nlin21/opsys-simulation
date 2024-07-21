#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

double next_exp(double lambda, int bound) {
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

	if (lambda <= 0) {
		fprintf(stderr, "ERROR: invalid lambda value\n");
		abort();
	}

	if (bound <= 0) {
		fprintf(stderr, "ERROR: invalid bound value\n");
		abort();
	}

	// 0 if CPU-bound, 1 if IO-bound
	bool type = 0;

	double CPU_total_CPU_burst_time = 0;
	double CPU_total_IO_burst_time = 0;
	int CPU_num_CPU_burst = 0;
	int CPU_num_IO_burst = 0;

	double IO_total_CPU_burst_time = 0;
	double IO_total_IO_burst_time = 0;
	int IO_num_CPU_burst = 0;
	int IO_num_IO_burst = 0;

	printf("<<< PROJECT PART I\n");
	if (n_CPU == 1) {
		printf("<<< -- process set (n=%d) with 1 CPU-bound process\n", n);
	} else {
		printf("<<< -- process set (n=%d) with %d CPU-bound processes\n", n, n_CPU);
	}
	printf("<<< -- seed=%d; lambda=%.6f; bound=%d\n", seed, lambda, bound);

	for (int i = 0; i < n; ++i) {
		if (i >= n_CPU) {
			type = 1;
		}

		int letter = 65 + (i / 10);
		int number = i % 10;

		// Step 1
		int arrival_time = floor(next_exp(lambda, bound));

		// Step 2
		double uniform_distribution = drand48();
		int CPU_bursts = ceil(uniform_distribution * 32);

		if (type == 0) {
			if (CPU_bursts == 1) {
				printf("CPU-bound process %c%i: arrival time %ims; 1 CPU burst:\n", letter, number, arrival_time);
			} else {
				printf("CPU-bound process %c%i: arrival time %ims; %i CPU bursts:\n", letter, number, arrival_time, CPU_bursts);
			}
		} else {
			if (CPU_bursts == 1) {
				printf("I/O-bound process %c%i: arrival time %ims; 1 CPU burst:\n", letter, number, arrival_time);
			} else {
				printf("I/O-bound process %c%i: arrival time %ims; %i CPU bursts:\n", letter, number, arrival_time, CPU_bursts);
			}
		}

		// Step 3
		for (int j = 0; j < CPU_bursts; ++j) {
			int CPU_burst_time = ceil(next_exp(lambda, bound));
			int IO_burst_time = -1;
			if (j < CPU_bursts-1) {	// do not generate IO burst time for the last CPU burst
				IO_burst_time = ceil(next_exp(lambda, bound));
			}
			if (type == 0) {	// CPU-bound
				CPU_burst_time *= 4;
				CPU_total_CPU_burst_time += CPU_burst_time;
				CPU_num_CPU_burst++;
				if (IO_burst_time != -1) {
					CPU_total_IO_burst_time += IO_burst_time;
					CPU_num_IO_burst++;
				}
			} else {			// IO-bound
				IO_burst_time *= 8;
				IO_total_CPU_burst_time += CPU_burst_time;
				IO_num_CPU_burst++;
				if (IO_burst_time != -8) {
					IO_total_IO_burst_time += IO_burst_time;
					IO_num_IO_burst++;
				}
			}

			if (j == CPU_bursts-1) {
				printf("==> CPU burst %ims\n", CPU_burst_time);
			} else {
				printf("==> CPU burst %ims ==> I/O burst %ims\n", CPU_burst_time, IO_burst_time);
			}
		}
	}

	double CPU_avg_CPU_burst_time = 0;
	if (CPU_num_CPU_burst != 0) {
		CPU_avg_CPU_burst_time = CPU_total_CPU_burst_time / CPU_num_CPU_burst;
	}

	double IO_avg_CPU_burst_time = 0;
	if (IO_num_CPU_burst != 0) {
		IO_avg_CPU_burst_time = IO_total_CPU_burst_time / IO_num_CPU_burst;
	}
	
	double avg_CPU_burst_time = 0;
	if (CPU_num_CPU_burst + IO_num_CPU_burst != 0) {
		avg_CPU_burst_time = (CPU_total_CPU_burst_time + IO_total_CPU_burst_time) / (CPU_num_CPU_burst + IO_num_CPU_burst);
	}

	double CPU_avg_IO_burst_time = 0;
	if (CPU_num_IO_burst != 0) {
		CPU_avg_IO_burst_time = CPU_total_IO_burst_time / CPU_num_IO_burst;
	}

	double IO_avg_IO_burst_time = 0;
	if (IO_num_IO_burst != 0) {
		IO_avg_IO_burst_time = IO_total_IO_burst_time / IO_num_IO_burst;
	}
	
	double avg_IO_burst_time = 0;
	if (CPU_num_IO_burst + IO_num_IO_burst != 0) {
		avg_IO_burst_time = (CPU_total_IO_burst_time + IO_total_IO_burst_time) / (CPU_num_IO_burst + IO_num_IO_burst);
	}

	FILE *file;
	file = fopen("simout.txt", "w");
	fprintf(file, "-- number of processes: %d\n", n);
	fprintf(file, "-- number of CPU-bound processes: %d\n", n_CPU);
	fprintf(file, "-- number of I/O-bound processes: %d\n", n - n_CPU);
	fprintf(file, "-- CPU-bound average CPU burst time: %.3f ms\n", ceil(CPU_avg_CPU_burst_time * 1000) / 1000);
	fprintf(file, "-- I/O-bound average CPU burst time: %.3f ms\n", ceil(IO_avg_CPU_burst_time * 1000) / 1000);
	fprintf(file, "-- overall average CPU burst time: %.3f ms\n", ceil(avg_CPU_burst_time * 1000) / 1000);
	fprintf(file, "-- CPU-bound average I/O burst time: %.3f ms\n", ceil(CPU_avg_IO_burst_time * 1000) / 1000);
	fprintf(file, "-- I/O-bound average I/O burst time: %.3f ms\n", ceil(IO_avg_IO_burst_time * 1000) / 1000);
	fprintf(file, "-- overall average I/O burst time: %.3f ms\n", ceil(avg_IO_burst_time * 1000) / 1000);
	fclose(file);

	return EXIT_SUCCESS;
}