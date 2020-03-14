#ifndef INPUT

#define INPUT
#include <string.h>
#include "tsp.h"
#define LINE "----------------------------------------------------------------------------------\n"
#define STAR_LINE "**********************************************************************************\n"
#define LINE_SIZE 180

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
	Deallocation of a tsp instance 
	\param tsp_in reference to tsp instance structure
*/
void dealloc_inst(tsp_instance*);

/**
	Print all command line parameters
*/
void help();

#endif 

