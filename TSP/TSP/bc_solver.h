/**
	@file bc_solver.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of CPLEX solutions with Branch & Cut.
*/

#ifndef BC_SOLVER
#define BC_SOLVER

#include "cplex_solver.h"

/**
	@brief Branch&Cut solver with lazy callbacks.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param tsp_in reference to tsp instance structure
	@param succ array of the successor of each node
	@param comp array with the component of each node
	@param general 0 if solved with directly lazy callbacks, 1 if solved with general callbacks
*/
void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general);

/**
	@brief Lazy Callback for adding of sec constraints.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param cbdata pointer at specific information for the callback
	@param wherefrom where the callback is called in the optimization
	@param cbhandle point to information by the user
	@param useraction_p specific what to do at the end of the callback
*/
static int CPXPUBLIC sec_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);

/**
	@brief Lazy Callback for adding of sec constraints (using general callbacks).
	@param context pointer to the context of the callback
	@param contextid specific the context in which call the callback
	@param cbhandle argument pass to the callback in the installation
*/
static int CPXPUBLIC sec_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);

/**
	@brief Function that adds constraints, one for each connected components.
	@param env pointer to the ENV structure, used by CPLEX solver
	@param tsp_in reference to tsp instance structure
	@param x_star solution in this node of the tree
	@param cbdata pointer at specific information for the callback
	@param wherefrom where the callback is called in the optimization
*/
int sec_bc_constraint(CPXENVptr env, tsp_instance* tsp_in, double* x_star, void* cbdata, int wherefrom);

/**
	@brief Function that adds constraints, one for each connected components.
	@param context context of the function
	@param tsp_in reference to tsp instance structure
	@param x_star solution in this node of the tree
*/
int sec_bc_constraint_general(CPXCALLBACKCONTEXTptr context, tsp_instance* tsp_in, double* x_star);

/**
	@brief Function that adds constraints for the hard fixing method.
	@param tsp_in reference to tsp instance structure
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param x_best solution in this node of the tree
	@param percentage how many edge to fix

*/
void cplex_change_coeff(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp, double* x_best, int percentage);

/**
	@brief Function that adds constraints for the soft fixing method.
	@param tsp_in reference to tsp instance structure
	@param env pointer to the ENV structure, used by CPLEX solver
	@param lp pointer to the LP structure, used by CPLEX solver
	@param x_best solution in this node of the tree
	@param freedom how many edge cplex can change in the next solution
*/
void local_branching(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp, double* x_best, int freedom);

static int CPXPUBLIC patching_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);

int patching(CPXENVptr env, tsp_instance* tsp_in, double* x_star, double objval, void* cbdata, int wherefrom);

static int CPXPUBLIC heuristic_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p);

static int CPXPUBLIC set_incumbent(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double objval, double* x, int* isfeas_p, int* useraction_p);

#endif
