#ifndef TSP

#define TSP
#include <string.h>

#define NUM_COMMANDS 8
#define DEADLINE_MAX 1000000

typedef struct
{
	//Input (Graph structure) 
	int num_nodes;
	double* x_coords;
	double* y_coords;

	//Parameters
	char input[500];
	long deadline;

} tsp_instance;

#endif
