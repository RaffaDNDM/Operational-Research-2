#ifndef HEURISTIC
#define HEURISTIC

#include "tsp.h"
#include "utility.h"
#include <time.h>
#include <float.h>
#include <math.h>

#define GRASP 
#define MAX_LOCAL_MINS 200
#define MAX_NUM_ITERATIONS 1000
#define FIRST_SOLUTION 1		// 0 = nearest neighborhood algorithm, 1 = insertion algorithm


/**
	Heuristic solver, call the select algorithm
	\param tsp_in reference to tsp instance structure
*/
void heuristic_solver(tsp_instance* tsp_in);

/**
	Insertion algorithm
	\param tsp_in reference to tsp instance structure
*/
void insertion(tsp_instance* tsp_in, int* visited_nodes, double* best_cost);

/**
	Nearest neighborhood algorithm
	\param tsp_in reference to tsp instance structure
*/
void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes, double* best_cost);

void min_cost(tsp_instance* tsp_in, int* nodes, int i, double* min_dist, int* best);

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost);

void greedy_refinement(tsp_instance* tsp_in, int* visited_nodes, double * cost);

void vns(tsp_instance* tsp_instance, int* visited_nodes, double* best_cost);

void update_solution(int* visited_nodes, double* sol, int num_nodes);

void succ_construction(int* visited_nodes, int* succ, int num_nodes);

void tabu_search(tsp_instance* tsp_in, int* visited_nodes, double* best_cost);

void add_element(int* list1, int* list2, int dimension, int element1, int element2);

void greedy_refinement_for_tabu_search(tsp_instance* tsp_in, int* visited_nodes, int** tabu_list, int tenure, int* num_tabu_edges, double* cost);
#endif
