#ifndef CPLEX_SOLVER

#define CPLEX_SOLVER

#include "utility.h"
#include <cplex.h>
#include <time.h>

//Set the use of lazy constraint
#define LAZY_CONSTRAINTS 1

//Check if the solution is an integer vector and respect constraaints
#define SOLUTION_CORRECTNESS 1

//Use metaheuristic pre-solution in LOOP
//#define METAHEURISTIC 1
//Approximation value for the cplex result values in pre-solution
#define EPS 1e-5

//Do all the combinations of CPLEX parameters
//#define ALL_PARAM_COMBINATIONS 1
#define NODE_LIMIT 3
#define EPS_GAP 0.1
#define	SOL_LIMIT 5 
#define SEED 500

//Number of ms that compose a sec (to convert second to ms)
#define TIME_SCALE 1000

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

/**
	MTZ definition of the model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void mtz_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/*
	Define the tour of the solution of MTZ solver
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void mtz_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp);

/**
	CPLEX definition of the model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void gg_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in);

/*
	Define the tour of the GG solution 
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void gg_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int compact_xpos(tsp_instance* tsp_in, int i, int j);

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
	Print state of an iteration of the loop algorithm
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param n_comp number of components in this iteration of the loop solver
	\param start time in which this iteration begins
	\param end time in which this iteration ends
*/
void print_state(CPXENVptr env, CPXLPptr lp, int ncomps, time_t start, time_t end);

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
