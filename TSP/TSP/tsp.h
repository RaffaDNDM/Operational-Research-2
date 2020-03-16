#ifndef TSP

#define TSP
#include <string.h>

#define NUM_COMMANDS 8
#define DEADLINE_MAX 1000000
#define VERBOSE 150
#define CAST_PRECISION 0.4999999999
#define NUM_ALGS 1

typedef struct
{
	//Input (Graph structure) 
	int num_nodes;
	double* x_coords;
	double* y_coords;

	//Parameters
	char input[500];
	long deadline;
	int alg;
	int integerDist;

	//Output
	int* sol;
	int bestCostI;
	double bestCostD;

} tsp_instance;

/**
	Compute the solution of the TSP problem 
	looking to the algorithm chosen by user on command line
	0 no permutation (default algorithm)
	\param tsp_in reference to tsp instance structure
*/
void solution(tsp_instance*);

/**
	Evaluate the cost of the solution already computed
	\param tsp_in reference to tsp instance structure
*/
void evaluate_sol(tsp_instance*);

/**
	Compute the distance between two nodes
	\param node1 index of first node (index of the node specified in TSP file)
	\param node2 index of second node (index of the node specified in TSP file)
	\param tsp_in reference to tsp instance structure
	\param dist container of the distance computed (coherent with costInt)
	\param costInt 1 if you want intenger distance, 0 otherwise  
*/
void dist(int, int, tsp_instance*, void*, int);

#endif

