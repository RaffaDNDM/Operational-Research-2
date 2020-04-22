#ifndef CPLEX_SOLVER

#define CPLEX_SOLVER

#include <cplex.h>
#include <time.h>
#include "utility.h"

//Number of ms that compose a sec (to convert second to ms)
#define TIME_SCALE 1000

//Set the use of lazy constraint
#define LAZY_CONSTRAINTS 1

//Check if the solution is an integer vector and respect constraaints
#define SOLUTION_CORRECTNESS 1

//Use metaheuristic pre-solution in LOOP
#define METAHEURISTIC 1
//Approximation value for the cplex result values in pre-solution
#define EPS 1e-5

//Do all the combinations of CPLEX parameters
#define ALL_PARAM_COMBINATIONS 1
#define NODE_LIMIT 3
#define EPS_GAP 0.1
#define	SOL_LIMIT 5 
#define SEED 500



/**
	CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void cplex_solver(tsp_instance* tsp_in);

/**
	CPLEX definition of the model
	\param tsp_in reference to tsp instance structure
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
*/
void cplex_build_model(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp);

/*
	Define the tour of the non-compact solution
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void cplex_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps);

/*
	Plot the CPLEX solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void cplex_plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int cplex_xpos(tsp_instance* tsp_in, int i, int j);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int compact_xpos(tsp_instance* tsp_in, int i, int j);

#endif
