#ifndef CPLEX_SOLVER

#define CPLEX_SOLVER

#include "utility.h"
#include <cplex.h>

#define SOLUTION_CORRECTNESS 1

#ifdef SOLUTION_CORRECTNESS
	
	//Approximation value for the cplex result values
	#define EPS 1e-5

#endif

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
	\param x array of the point in the solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void cplex_define_tour(double*, tsp_instance*, int*, int*, int*);

/*
	Plot the CPLEX solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void plot_cplex(tsp_instance*, int*, int*, int*);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int xpos_cplex(tsp_instance*, int, int);

/**
	CPLEX definition of the model
	\param tsp_in reference to tsp instance structure
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
*/
void mtz_build_model(tsp_instance*, CPXENVptr, CPXLPptr);


/**
	CPLEX definition of the model
	\param tsp_in reference to tsp instance structure
	\param env pointer to the ENV structure, used by CPLEX solver
	\param lp pointer to the LP structure, used by CPLEX solver
*/
void gg_build_model(tsp_instance*, CPXENVptr, CPXLPptr);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int compact_xpos(tsp_instance*, int, int);

/*
	Define the tour fo the solution
	\param x array of the point in the solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void compact_define_tour(double*, tsp_instance*, int*, int*);


#endif