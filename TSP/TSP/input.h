#ifndef INPUT

#define INPUT

#include "utility.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

/**
	Parser of command line
	\param argv array of command line parameters
	\param argc number of command line parameters
	\param tsp_in reference to tsp instance structure 
*/
void parse_cmd(char** argv, int argc, tsp_instance* tsp_in);

/**
	Parser of tsp file
	\param tsp_in reference to tsp instance structure 
*/
void parse_file(tsp_instance* tsp_in);

/**
	Parser of tsp file
	\param tsp_in reference to tsp instance structure
	\param alg_string string that defines types of algorithms
	\param in_main 1 if in cycle, 0 in parser  of command line
*/
void select_alg(tsp_instance* tsp_in, char* alg_string, int in_main);

/**
	Deallocaion of a tsp instance 
	\param tsp_in reference to tsp instance structure
*/
void dealloc_inst(tsp_instance* tsp_in);

/**
	Print all command line parameters
*/
void help();

#endif 

