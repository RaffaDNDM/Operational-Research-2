#ifndef HEURISTIC
#define HEURISTIC

#include "tsp.h"
#include "utility.h"
#include <time.h>
#include <float.h>

/**
	Heuristic solver, call the select algorithm
	\param tsp_in reference to tsp instance structure
*/
void heuristic_solver(tsp_instance* tsp_in);

/**
	Insertion algorithm
	\param tsp_in reference to tsp instance structure
*/
void insertion(tsp_instance* tsp_in);

/**
	Nearest neighborhood algorithm
	\param tsp_in reference to tsp instance structure
*/
void nearest_neighborhood(tsp_instance* tsp_in);


#endif
