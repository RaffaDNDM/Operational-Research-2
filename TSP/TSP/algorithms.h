#include "tsp.h"
#include <cplex.h>

#define LP_FILENAME "model.lp"
#define NAME_SIZE 30

/**
	Default algorithm for TSP
	\param tsp_in reference to tsp instance structure
*/
void default_alg(tsp_instance*);

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


void path_tsp(double* x, int size, tsp_instance* tsp_in);