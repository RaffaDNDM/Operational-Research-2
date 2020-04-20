#ifndef CPLEX_SOLVER

#define CPLEX_SOLVER

#include "utility.h"
#include <cplex.h>
#include <time.h>

#define LAZY_CONSTRAINTS 1
#define SOLUTION_CORRECTNESS 1

//#define METAHEURISTIC 1
//#define ALL_PARAM_COMBINATIONS 1
#define NODE_LIMIT 3
#define EPS_GAP 0.1
#define	ORDER_SOL 5 
#define SEED 500
#define TIME_SCALE 1000

//Approximation value for the cplex result values
#define EPS 1e-5

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
	Define the tour fo the solution
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void cplex_define_tour(tsp_instance*, double*, int*, int*, int*);

/*
	Plot the CPLEX solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void cplex_plot(tsp_instance*, int*, int*, int*);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int cplex_xpos(tsp_instance*, int, int);

/**
	CPLEX definition of the model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void mtz_build_model(CPXENVptr, CPXLPptr, tsp_instance*);

/*
	Define the tour fo the solution
		\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void mtz_define_tour(tsp_instance*, double*, int*, int*);

/**
	CPLEX definition of the model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
*/
void gg_build_model(CPXENVptr, CPXLPptr, tsp_instance*);

/*
	Define the tour fo the solution
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void gg_define_tour(tsp_instance*, double*, int*, int*);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int compact_xpos(tsp_instance*, int, int);

/**
	Benders definition of the model and resolution
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param succ
	\param comp
*/
void loop_solver(CPXENVptr, CPXLPptr, tsp_instance*, int*, int *);

/**
	CPLEX definition of the model
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
	\param tsp_in reference to tsp instance structure
	\param succ
	\param comp
*/
void add_sec_constraint(CPXENVptr, CPXLPptr, tsp_instance*, int*, int);


void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general);

void print_state(CPXENVptr, CPXLPptr, int, time_t, time_t);

static int CPXPUBLIC sec_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p);
static int CPXPUBLIC sec_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle);

int sec_bc_constraint(CPXENVptr env, tsp_instance* tsp_in, double* x_star, void* cbdata, int wherefrom);
int sec_bc_constraint_general(CPXCALLBACKCONTEXTptr context, tsp_instance* tsp_in, double* x_star);
#endif
