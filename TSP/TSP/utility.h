/**
	@file utility.h
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Header of utility functions.
*/

#ifndef UTILITY
#define UTILITY

#include "tsp.h"

#define LINE "----------------------------------------------------------------------------------\n"
#define STAR_LINE "**********************************************************************************\n"

#define LINE_SIZE 180 //size of the line gets from the input
#define NAME_SIZE 30 //size of the variable name in the cplex model


//Check if the solution is an integer vector and respect constraints
#define SOLUTION_CORRECTNESS 1

//Approximation value for the cplex result values in pre-solution
#define EPS 1e-5

#define GP_DEFAULT_STYLE "style.txt" //name of the style file needed by Gnuplot
#define DEFAULT_DAT "solution.dat" //name of the file where to print the solution of the default algoritm
#define CPLEX_DAT "solutionCPLEX.dat" //name of the file where to print the cplex solution
#define GP_CPLEX_STYLE "styleCPLEX.txt" //name of the style file for the cplex solution needed by Gnuplot
#define LP_FILENAME "model.lp" //name fo the file that conteins the cplex model
#define GNUPLOT_EXE "%GNUPLOT%/bin/gnuplot.exe -persistent" //path of  where is gnuplot.exe
#define PLOT_HEURISTIC_DAT "plot_cost.dat"
#define HEURISTIC_STYLE "styleHEURISTIC.txt"

#define CAST_PRECISION 0.4999999999 //quantity to add to correctly approximate double into int

#define DEFAULT_SOLLIM_VALUE 2147483647

//Colors
#define BLUE "\033[1;34m"
#define GREEN "\033[1;32m"
#define RED "\033[1;31m"
#define WHITE "\033[0m"
#define YELLOW "\033[1;33m"
#define CYAN "\033[1;36m"

/**
	@brief Compute the distance between two nodes, looking to the specified way of computing distances in tsp_in.
	@param node1 index of first node (index of the node specified in TSP file)
	@param node2 index of second node (index of the node specified in TSP file)
	@param tsp_in reference to tsp instance structure
	@param dist container of the distance computed (coherent with costInt)
*/
void dist(int, int, tsp_instance*, void*);

/**
	@brief Print value of the solution of TSP problem.
	@param tsp_in reference to tsp instance structure
*/
void print_cost(tsp_instance* tsp_in);

/**
	@brief Return the position of the element (i,j) in the matrix of corresponding edges, created from the nodes in the graph.
	@param tsp_in reference to tsp instance structure
	@param i first index
	@param j second index
*/
int xpos(tsp_instance* tsp_in, int i, int j);

int generic_xpos(int i, int j, int num_nodes);

/**
	@brief Define the tour of the non-compact solution.
	@param tsp_in reference to tsp instance structure
	@param x array of the point in the solution
	@param succ array of the successor of each node
	@param comp array with the component of each node
	@param n_comps number of components
*/
void define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps);

/**
	@brief Plot the solution.
	@param tsp_in reference to tsp instance structure
	@param succ array of the successor of each node
	@param comp array with the component of each node
	@param n_comps number of components
*/
void plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps);

void cost_plot_definition(tsp_instance* tsp_in);


#endif
