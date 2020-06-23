/**
	@file heuristic.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of Heuristic solvers.
*/

#ifndef HEURISTIC
#define HEURISTIC

#include "tsp.h"
#include "utility.h"
#include <time.h>
#include <float.h>
#include <math.h>

#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>
//#define MULTI_START //Comment or not if you want multistart or not

#define STEP_SEED 100
#define NUM_MULTI_START 12 //Number of threads in multi start and for construction in Genetic
#define POPULATION_SIZE 12000 //Size of the population
#define NUM_WORST_MEMBERS 1000 //Size of the buffer with indices of worst members in the population
#define FIXED_TIME_MS
//#define UNIFORM_PROB

typedef struct
{
	tsp_instance* tsp_in; //pointer to tsp instance
	int* succ; //pointer to the best visited nodes sequence
	int seed; //seed used by each thread
	time_t start;
}thread_args;

typedef struct
{
	tsp_instance* tsp_in; //pointer to tsp instance
	int** members; //members of the population
	int num_members;
	int first_index;
	double* fitnesses; //costs of each
	int* num_instances;
	double* sum_fitnesses;
	double* sum_prob;
	int* best_index;
}construction_args;

typedef struct
{
	int end_list;//-1
	int start_list; //0
} tabu_list_params;

#define GRASP
#define MAX_LOCAL_MINS 200
#define MAX_NUM_ITERATIONS 1000
#define CONSTRUCTION_TYPE 0	// 0 = nearest neighborhood algorithm, 1 = insertion algorithm
#define REACTIVE //define for use the reactive tabu search
#define MAX_NUM_EPOCHS 100

/**
	@brief Compute the solution of the instance invoking the correct function.
	@param param pointer to a structure with the needed information for the computation
*/

void* computeSolution(void* param);

/**
	@brief Heuristic solver, call the select algorithm.
	@param tsp_in reference to tsp instance structure
*/
void heuristic_solver(tsp_instance* tsp_in);

/**
	@brief Insertion algorithm.
	@param tsp_in reference to tsp instance structure
*/
void insertion(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, int seed, int first_node);

/**
	@brief Nearest neighborhood algorithm.
	@param tsp_in reference to tsp instance structure
*/
void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, int seed, int first_node);

/**
	@brief find the free node at the min distance.
	@param tsp_in reference to tsp instance structure
	@param nodes array with 1 if the node is already visited, 0 otherwise 
	@param i starting node
	@param min_dist at the end, it will contain the minimum distance find 
	@param best at the end, it will contain the position of the node find
	@param seed random seed for the GRASP option
*/

void min_cost(tsp_instance* tsp_in, int* nodes, int i, double* min_dist, int* best, int seed);

/**
	@brief Compute the minimum extra-mileage.
	@param tsp_in reference to tsp instance structure
	@param count number of nodes already in the solution
	@param visited_nodes array with sequence of the visited nodes
	@param node1 array with the first node of the edges in the tour
	@param node2 array with the second node of the edges in the tour
	@param costs array with the cost of each edge in the tour
	@param i_best node in the tour selected
	@param k_best node to insert
	@param best_cost_h variation of the cost 
	@param best_cost final cost of the tour
	@param seed random seed for the GRASP option
*/

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost, int seed);

/**
	@brief Compute a 2-opt refinement to the actual solution
	@param tsp_in reference to tsp instance structure
	@param visited_nodes array with sequence of the visited nodes
	@param cost cost of the solutoion
*/

void greedy_refinement(tsp_instance* tsp_in, int* visited_nodes, double * cost);

//void vns(tsp_instance* tsp_instance, int* visited_nodes, double* best_cost, double deadline);

/**
	@brief Compute the hybrid VNS algorithm
	@param tsp_in reference to tsp instance structure
	@param visited_nodes array with sequence of the visited nodes
	@param best_cost cost of the solution
	@param deadline time limit
*/

void hybrid_vns(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, double deadline);

/**
	@brief Find the k-opt sequence at the minimun distance
	@param tsp_in reference to tsp instance structure
	@param kopt_visited_nodes sequence of the visited nodes of the k-opt sequence found
	@param kopt_cost cost of the solution found
	@param k distance between the two nodes to swap
	@param inverse_costs array with the inverse costs
	@param inverse_costs_sum sum of the inverse costs
*/

int min_kopt_sequence(tsp_instance* tsp_in, int* kopt_visited_nodes, double* kopt_cost, int k, double** inverse_costs, double* inverse_costs_sum);

/**
	@brief Compute a new random solution
	@param tsp_in reference to tsp instance structure
	@param local_min_visited_nodes sequence of the visited nodes of the local minimun
	@param local_cost cost of the local minimun
	@param inverse_costs array with the inverse costs
	@param inverse_costs_sum sum of the inverse costs
*/

int new_random_sol(tsp_instance* tsp_in, int* local_min_visited_nodes, double* local_cost, double** inverse_costs, double* inverse_costs_sum);

/**
	@brief Update the array with the solution.
	@param visited_nodes array with sequence of the visited nodes
	@param sol array with the solution
	@param num_nodes number of nodes in the problem
*/

void update_solution(int* visited_nodes, double* sol, int num_nodes);

/**
	@brief Build the vector of the successors
	@param visited_nodes array with sequence of the visited nodes
	@param succ array of the successors
	@param num_nodes number of nodes in the problem
*/

void succ_construction(int* visited_nodes, int* succ, int num_nodes);

/**
	@brief Compute the Tabu Search algorithm
	@param tsp_in reference to tsp instance structure
	@param visited_nodes array with sequence of the visited nodes
	@param best_cost cost of the solution
	@param deadline time limit
*/

void tabu_search(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, double deadline);

/**
	@brief Add an element to the tabu list
	@param list1 array with the firsts nodes of the edges in the tabu list
	@param list2 array with the seconds nodes of the edges in the tabu list
	@param dimension max dimension of the tabu list
	@param element1 first node of the edge to add 
	@param element2 second node of the edge to add
	@param with_reduction 1 if the tabu list has to be reduced, 0 otherwise
	@param logically_full 1 if the tabu list cannot be enlarged, 0 otherwise
	@param param information needed to manage the list
*/

void add_element(int* list1, int* list2, int dimension, int element1, int element2, int with_reduction, int logically_full, tabu_list_params* param);

/**
	@brief 2-opt refinement respecting a tabu list.
	@param tsp_in reference to tsp instance structure
	@param tabu_list list of the tabu edges
	@param param information needed to manage the list
	@param max_tenure max dimension of the tabu list
	@param min_tenure min dimension of the tabu list
	@param num_tabu_edges number of edges in the tabu list
	@param cost cost of the solution
*/

void greedy_refinement_for_tabu_search(tsp_instance* tsp_in, int* visited_nodes, int** tabu_list, tabu_list_params* param,int max_tenure, 
	int min_tenure, int* num_tabu_edges, double* cost);

/**
	@brief Compute the Simulated Annealing algorithm.
	@param tsp_in reference to tsp instance structure
	@param visited_nodes array with sequence of the visited nodes
	@param best_cost cost of the solution
	@param deadline time limit
*/

void simulated_annealing(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, double deadline);

/**
	@brief Manage the genetic algorithm.
	@param tsp_in reference to tsp instance structure
*/

void genetic_solver(tsp_instance* tsp_in);

/** 
	@brief Build the first element of the population fo the genetic algorithm.
	@param param information needed for the construction
*/

void construction(void* param);

/**
	@brief Compute the evolution of the genetic algorithm.
	@param tsp_in reference to tsp instance structure
	@param members visited nodes if each member of the population
	@param fitnesses cost of the solution of each member of the population
	@param best_index index of the members with the minimum cost
	@param worst_members array with the members with the maximum cost
	@param sum_prob sum of the probabilities
	@param sum_fitnesses sum of the fitnesses 
	@param start starting time of the computation
*/

void evolution(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index, int* worst_members, double* sum_prob, double* sum_fitnesses, time_t start);

/** 
	@brief Compute the crossover in the evolution of the genetic algorithm.
	@param tsp_in reference to tsp instance structure
	@param members visited nodes if each member of the population
	@param fitnesses cost of the solution of each member of the population
	@param best_index index of the members with the minimum cost
	@param worst_members array with the members with the maximum cost
	@param sum_prob sum of the probabilities
	@param sum_fitnesses sum of the fitnesses
	@param seed random seed
	@param index first_index of worst_members
*/

void crossover(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index, int* worst_members, double* sum_prob, double* sum_fitnesses, int seed, int* index);

/**
	@brief Compute the mutation in the evolution of the genetic algorithm.
	@param tsp_in reference to tsp instance structure
	@param members visited nodes if each member of the population
	@param fitnesses cost of the solution of each member of the population
	@param best_index index of the members with the minimum cost
	@param worst_members array with the members with the maximum cost
	@param sum_prob sum of the probabilities
	@param sum_fitnesses sum of the fitnesses
	@param seed random seed
	@param index  first_index of worst_members
*/

void mutation(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index, int* worst_members, double* sum_prob, double* sum_fitnesses, int seed, int* index);

/**
	@brief Update the array with the worst members of the population.
	@param fitnesses cost of the solution of each member of the population
	@param worst_members array with the members with the maximum cost
*/

void update_worst(double* fitnesses, int* worst_members);

/**
	@brief Check if an edge is in the tabu list, return 1 if the edge is forbidden, 0 otherwise
	@param tabu_list  list of the tabu edges
	@param tenure dimension of the tabu list
	@param node1 first node of the edge to check
	@param node2 second node of the edge to check
*/

int check_tabu_list(int** tabu_list, int tenure, int node1, int node2);

/**
	@brief Compute a 2-opt move respecting a tabu list, return the minimum variation of the cost
	@param tsp_in reference to tsp instance structure
	@param succ array with the successors
	@param tabu_list list of the tabu edges
	@param tenure dimension of the tabu list
	@param param information needed to manage the list
*/

double move2opt_for_tabu_search(tsp_instance* tsp_in, int* succ, int** tabu_list, int* tenure, tabu_list_params* params);


#endif

