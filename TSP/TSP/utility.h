#ifndef UTILITY

#define UTILITY

#include "tsp.h"

#define LINE "----------------------------------------------------------------------------------\n"
#define STAR_LINE "**********************************************************************************\n"

#define LINE_SIZE 180 //size of the line gets from the input
#define NAME_SIZE 30 //size of the variable name in the cplex model

#define GP_DEFAULT_STYLE "style.txt" //name of the style file needed by Gnuplot
#define DEFAULT_DAT "solution.dat" //name of the file where to print the solution of the default algoritm
#define CPLEX_DAT "solutionCPLEX.dat" //name of the file where to print the cplex solution
#define GP_CPLEX_STYLE "styleCPLEX.txt" //name of the style file for the cplex solution needed by Gnuplot
#define LP_FILENAME "model.lp" //name fo the file that conteins the cplex model 
#define GNUPLOT_EXE "%GNUPLOT%/bin/gnuplot.exe -persistent" //path of  where is gnuplot.exe 

#define CAST_PRECISION 0.4999999999 //quantity to add to correctly approximate double into int

#define DEFAULT_SOLLIM_VALUE 2147483647 

/**
	Compute the distance between two nodes
	(looking to the specified way of computing distances in *tsp_in)
	\param node1 index of first node (index of the node specified in TSP file)
	\param node2 index of second node (index of the node specified in TSP file)
	\param tsp_in reference to tsp instance structure
	\param dist container of the distance computed (coherent with costInt)
*/
void dist(int, int, tsp_instance*, void*);

/*
	Print value of the solution of TSP problem
	\param tsp_in reference to tsp instance structure
*/
void print_cost(tsp_instance* tsp_in);

/*
	Return the position of the element (i,j) in the matrix
	of corresponding edges, created from the nodes in the graph
	\param tsp_in reference to tsp instance structure
	\param i first index
	\param j second index
*/
int xpos(tsp_instance* tsp_in, int i, int j);

/*
	Define the tour of the non-compact solution
	\param tsp_in reference to tsp instance structure
	\param x array of the point in the solution
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps);

/*
	Plot the solution
	\param tsp_in reference to tsp instance structure
	\param succ array of the successor of each node
	\param comp array with the component of each node
	\param n_comps number of components
*/
void plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps);

#endif