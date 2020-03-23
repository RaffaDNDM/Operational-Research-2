#include "tsp.h"
#include <cplex.h>
#include "input.h"

#define NAME_SIZE 30

/**
	Default algorithm for TSP
	\param tsp_in reference to tsp instance structure
*/
void default_alg(tsp_instance*);

/**
	Compute the solution of the TSP problem
	looking to the algorithm chosen by user on command line
	1 no permutation (default algorithm)
	\param tsp_in reference to tsp instance structure
*/
void solution(tsp_instance*);

/**
	Evaluate the cost of the solution already computed
	(on default the distance is computed through double operations
	if user types -int or -i, it'computed through integer operations)
	\param tsp_in reference to tsp instance structure
*/
void evaluate_sol(tsp_instance*);

/**
	CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void cplex_solver(tsp_instance*);

/**
	CPLEX definition of the model
	\param tsp_in reference to tsp instance structure
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
*/
void cplex_build_model(tsp_instance*, CPXENVptr, CPXLPptr);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int xpos(tsp_instance*, int, int);

/*
	Print value of the solution of TSP problem
	\param tsp_in reference to tsp instance structure
*/
void print_sol(tsp_instance* tsp_in);

/**
	Plot solutions in 2D plane
	\param tsp_in reference to tsp instance structure
*/
void plot_solution(tsp_instance*);


void plot_cplex(double* x, int size, tsp_instance* tsp_in);