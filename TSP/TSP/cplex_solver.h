/**
	@file cplex_solver.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of CPLEX solvers.
*/

#ifndef CPLEX_SOLVER

#define CPLEX_SOLVER

#include <cplex.h>
#include <time.h>
#include "utility.h"

//Number of ms that compose a sec (to convert second to ms)
#define TIME_SCALE 1000

//Set the use of lazy constraint
//#define LAZY_CONSTRAINTS 1


//Use metaheuristic pre-solution in LOOP
//#define METAHEURISTIC 1


//Do all the combinations of CPLEX parameters
//#define ALL_PARAM_COMBINATIONS 1
#define NODE_LIMIT 3
#define EPS_GAP 0.1
#define	SOL_LIMIT 5
#define SEED 500



/**
	@brief CPLEX solver, call the select algorithm.
	@param tsp_in reference to tsp instance structure
*/
void cplex_solver(tsp_instance* tsp_in);

/**
	@brief CPLEX definition of the model.
	@param tsp_in reference to tsp instance structure
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
*/
void cplex_build_model(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp);

/**
	@brief Return the position of the element (i,j) in the matrix of corresponding edges, created from the nodes in the graph, in a compact model.
	@param tsp_in reference to tsp instance structure
	@param i first index
	@param j second index
*/
int compact_xpos(tsp_instance* tsp_in, int i, int j);


#endif
