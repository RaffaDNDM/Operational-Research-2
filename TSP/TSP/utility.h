#include "tsp.h"

#define LINE "----------------------------------------------------------------------------------\n"
#define STAR_LINE "**********************************************************************************\n"
#define LINE_SIZE 180

#define NAME_SIZE 30

#define GP_DEFAULT_STYLE "style.txt"
#define DEFAULT_DAT "solution.dat"
#define CPLEX_DAT "solutionCPLEX.dat"
#define GP_CPLEX_STYLE "styleCPLEX.txt"
#define LP_FILENAME "model.lp"
#define GNUPLOT_EXE "%GNUPLOT%/bin/gnuplot.exe -persistent"

#define CAST_PRECISION 0.4999999999

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