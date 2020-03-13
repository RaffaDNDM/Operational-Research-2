#ifndef INPUT

#define INPUT
#include <string.h>
#include "tsp.h"
#define LINE "-----------------------------------------------------\n"
#define LINE_SIZE 180

void parse_cmd(char**, int, tsp_instance*);
void help();
void parse_file(tsp_instance*);

#endif 

