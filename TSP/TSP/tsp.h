#ifndef TSP

#define TSP

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define DEADLINE_MAX 1000000
#define NUM_ALGS 2
#define NUM_COMMANDS 8
#define VERBOSE 150

typedef struct
{
	//Input (Graph structure) 
	int num_nodes;
	double* x_coords;
	double* y_coords;

	//Parameters
	char input[500];
	double deadline;
	int alg;
	int integerDist;
	int plot;
	int verbose;

	//Output
	int* sol;
	int bestCostI;
	double bestCostD;

} tsp_instance;

/**
	Compute the solution of the TSP problem
	looking to the algorithm chosen by user on command line
	1 no permutation (default algorithm)
	\param tsp_in reference to tsp instance structure
*/
void solution(tsp_instance*);

#endif

