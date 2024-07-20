#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>

#define ROWS 26
#define COLS 10

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

		if (p == 0) break;
	}

	if (p == 0) {	// CHILD
		printf("CHILD: Assigning process id %d to ID table\n", getpid());
		assign(getpid(), ID_TABLE);
	}

	return EXIT_SUCCESS;
}