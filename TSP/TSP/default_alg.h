#include "tsp.h"
#include "utility.h"

/**
	Default algorithm for TSP
	\param tsp_in reference to tsp instance structure
*/
void default_alg(tsp_instance*);

/**
	Evaluate the cost of the solution already computed
	(on default the distance is computed through double operations
	if user types -int or -i, it'computed through integer operations)
	\param tsp_in reference to tsp instance structure
*/
void evaluate_sol(tsp_instance*);

/**
	Plot solutions in 2D plane
	\param tsp_in reference to tsp instance structure
*/
void plot_solution(tsp_instance*);