#ifndef TSP

#define TSP

#include <string.h>
#include <assert.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define GP_DEFAULT_STYLE "style.txt"
#define DEFAULT_DAT "solution.dat"
#define CPLEX_DAT "solutionCPLEX.dat"
#define GP_CPLEX_STYLE "styleCPLEX.txt"
#define LP_FILENAME "model.lp"
#define GNUPLOT_EXE "%GNUPLOT%/bin/gnuplot.exe -persistent"

#define NUM_COMMANDS 8
#define DEADLINE_MAX 1000000
#define VERBOSE 150
#define CAST_PRECISION 0.4999999999
#define NUM_ALGS 2

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
	Compute the distance between two nodes
	(looking to the specified way of computing distances in *tsp_in)
	\param node1 index of first node (index of the node specified in TSP file)
	\param node2 index of second node (index of the node specified in TSP file)
	\param tsp_in reference to tsp instance structure
	\param dist container of the distance computed (coherent with costInt)
*/
void dist(int, int, tsp_instance*, void*);

#endif

