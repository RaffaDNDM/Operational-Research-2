#ifndef BC_SOLVER
#define BC_SOLVER

#include "cplex_solver.h"

/**
	Branch&Cut solver with lazy callbacks
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param general 0 if solved with directly lazy callbacks, 1 if solved with general callbacks
*/
void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general);

/**
	Lazy Callback for adding of sec constraints
	\param env pointer to the ENV structure, used by CPLEX solver
	\param cbdata
	\param wherefrom
	\param cbhandle
	\param useraction_p
*/
static int CPXPUBLIC sec_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);

/**
	Lazy Callback for adding of sec constraints (using general callbacks)
	\param context
	\param contextid
	\param cbhandle
*/
static int CPXPUBLIC sec_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);

/**
	Function that adds constraints, one for each connected components
	\param env pointer to the ENV structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param x_star solution in this node of the tree
	\param cbdata
	\param wherefrom
*/
int sec_bc_constraint(CPXENVptr env, tsp_instance* tsp_in, double* x_star, void* cbdata, int wherefrom);

/**
	Function that adds constraints, one for each connected components
	\context
	\param tsp_in reference to tsp instance structure
	\x_star solution in this node of the tree
*/
int sec_bc_constraint_general(CPXCALLBACKCONTEXTptr context, tsp_instance* tsp_in, double* x_star);

#endif