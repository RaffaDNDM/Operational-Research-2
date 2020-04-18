#ifndef TSP

#define TSP

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

//Maximum time for the computation of the solution (15 min)
#define DEADLINE_MAX 900 

//Number of algoritms implemented
#define NUM_ALGS 4

//#define NUM_COMMANDS 8

//default value of verbose
#define VERBOSE 40

#define DIR_LIST_PY "python3 list_dir.py"

#define ALG1 "Loop"
#define ALG2 "Branch&Bound"
#define ALG3 "MTZ"
#define ALG4 "GG"

#define PERF_PROF_ON 1
#define PERF_PROF_PY "python perfprof.py -D , -T 3600 -S 2 -M 20 perf_data.csv pp.pdf -P \"all instances, shift 2 sec.s\""

typedef struct
{
	//Input (Graph structure) 
	int num_nodes;
	double* x_coords;
	double* y_coords;
	int num_cols;

	//Parameters
	char input[500];
	double deadline;
	int alg;
	int which_alg[NUM_ALGS];
	int integerDist;
	int plot;
	int verbose;
	int node_lim;
	int sol_lim;
	double eps_gap;
	int seed;
	char dir[500];
	

	//Output
	int* sol;
	int bestCostI;
	double bestCostD;
	double execution_time;

} tsp_instance;

/**
	Compute the solution of the TSP problem
	looking to the algorithm chosen by user on command line
	1 no permutation (default algorithm)
	\param tsp_in reference to tsp instance structure
*/
void solution(tsp_instance*);

void set_params_and_solve(tsp_instance*);
#endif

