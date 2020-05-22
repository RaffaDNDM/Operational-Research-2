#ifndef HEURISTIC
#define HEURISTIC

#include "tsp.h"
#include "utility.h"
#include <time.h>
#include <float.h>

#define GRASP 1

/**
	Heuristic solver, call the select algorithm
	\param tsp_in reference to tsp instance structure
*/
void heuristic_solver(tsp_instance* tsp_in);

/**
	Insertion algorithm
	\param tsp_in reference to tsp instance structure
*/
void insertion(tsp_instance* tsp_in, int* visited_nodes);

/**
	Nearest neighborhood algorithm
	\param tsp_in reference to tsp instance structure
*/
void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes);

void min_cost(tsp_instance* tsp_in, int* nodes, int i, double* min_dist, int* best);

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost);

void greedy_refinement(tsp_instance* tsp_in, int* visited_nodes);

void vns(tsp_instance* tsp_in, int* visited_nodes);

void update_solution(tsp_instance* tsp_in, int* visited_nodes);
#endif
