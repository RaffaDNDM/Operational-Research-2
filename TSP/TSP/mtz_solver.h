/**
	@file mtz_solver.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of MTZ solver.
*/

#ifndef MTZ_SOLVER
#define MTZ_SOLVER

#include "cplex_solver.h"

void mtz_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/**
	@brief MTZ definition of the model.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param tsp_in reference to tsp instance structure
*/
void mtz_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/**
	@brief Define the tour of the solution of MTZ solver.
	@param tsp_in reference to tsp instance structure
	@param x array of the point in the solution
	@param succ array of the successor of each node
	@param comp array with the component of each node
*/
void mtz_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp);

#endif
