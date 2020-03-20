#ifndef INPUT

#define INPUT
#define CIAO 1

#include "tsp.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#define LINE "----------------------------------------------------------------------------------\n"
#define STAR_LINE "**********************************************************************************\n"
#define LINE_SIZE 180
//#define GNUPLOT_EXE "D:/Programs/Gnuplot/bin/gnuplot.exe -persistent"
#define GNUPLOT_EXE "D:/Programs/Gnuplot/bin/gnuplot.exe -persistent"
#define GNUPLOT_STYLE "style.txt"
#define SOLUTION_FILENAME "solution.dat"

/**
	Parser of command line
	\param argv array of command line parameters
	\param argc number of command line parameters
	\param tsp_in reference to tsp instance structure 
*/
void parse_cmd(char**, int, tsp_instance*);

/**
	Parser of tsp file
	\param tsp_in reference to tsp instance structure 
*/
void parse_file(tsp_instance*);

/**
	Plot solutions in 2D plane
	\param tsp_in reference to tsp instance structure
*/
void plot_solution(tsp_instance*);

/**
	Deallocation of a tsp instance 
	\param tsp_in reference to tsp instance structure
*/
void dealloc_inst(tsp_instance*);

/**
	Print all command line parameters
*/
void help();

#endif 

