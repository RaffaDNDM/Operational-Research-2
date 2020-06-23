/**
	@file gg_solver.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of Gavish Graves solver.
*/

#ifndef GG_SOLVER
#define GG_SOLVER

#include "cplex_solver.h"

/**
	@brief CPLEX solver with GG model.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param tsp_in reference to tsp instance structure
*/
void gg_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/**
	@brief CPLEX definition of the model.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param tsp_in reference to tsp instance structure
*/
void gg_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/**
	@brief Define the tour of the GG solution.
	@param tsp_in reference to tsp instance structure
	@param x array of the point in the solution
	@param succ array of the successor of each node
	@param comp array with the component of each node
	@param n_comps number of components
*/
void gg_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp);

#endif