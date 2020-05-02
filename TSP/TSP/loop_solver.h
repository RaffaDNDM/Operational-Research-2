#ifndef LOOP_SOLVER
#define LOOP_SOLVER

#include "cplex_solver.h"


/**
	Benders definition of the model and resolution
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	*/
void loop_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp);

/**
	Add all the sec constraints to the initial model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param comp array with the component of each node
	\param n_comps number of connected components in the solution
*/
void add_sec_constraint(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* comp, int n_comps);

/**
	Print state of an iteration of the loop algorithm
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param ncomps number of components in this iteration of the loop solver
	\param start time in which this iteration begins
	\param end time in which this iteration ends
*/
void print_state(CPXENVptr env, CPXLPptr lp, int ncomps, time_t start, time_t end);


#endif