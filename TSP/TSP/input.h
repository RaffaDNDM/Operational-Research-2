#ifndef INPUT

#define INPUT

#include "tsp.h"
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

