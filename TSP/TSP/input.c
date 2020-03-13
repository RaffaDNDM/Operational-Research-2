#include "tsp.h"
#include "input.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

void parse_cmd(char** argv, int argc, tsp_instance* tsp_in)
{
	//at leat one param + program name
	assert(argc < 2);
	char** commands = malloc(sizeof(char*)*NUM_COMMANDS);

	//verbose
	tsp_in->num_nodes = -1;
	tsp_in->num_nodes = DEADLINE_MAX;

	int def_deadline = 0;

	strcpy(tsp_in->input, "NULL");

	int i = 1;
	for (; i < argc; i++)
	{
		/*
		if(strncmp(argv[i], "-alg")==0)
			argv[++i]
		*/
		if ((strncmp(argv[i], "-d",2) == 0 || strncmp(argv[i], "-dead",5) == 0 || strncmp(argv[i], "-deadline",9) == 0) && !def_deadline)
		{
			def_deadline = 1;
			tsp_in->deadline = atol(argv[++i]);
			continue;
		}

		if ((strncmp(argv[i], "-f", 2) == 0 || strncmp(argv[i], "-file", 5) == 0 || strncmp(argv[i], "-input", 6) == 0))
		{
			if (strncmp(tsp_in->input, "NULL", 2))
			{
				strcpy(tsp_in->input, argv[++i]);
				continue;
			}
		}

		if ((strncmp(argv[i], "-help", 5) == 0))
		{
			help();
			continue;
		}
	}

	//VERBOSE
}

void help()
{
	printf("Help\n");
	printf(LINE);
	printf("Insert the file in input\n");

	printf("Insert the max time of the execution\n");

	printf("Ask for help about commands\n");

}


void parse_file(tsp_instance* tsp_in)
{
	FILE* f = fopen(tsp_in->input, "r");

	assert(f == NULL);
	
	tsp_in->num_nodes = -1;

	char line[LINE_SIZE];
	char* token;
	int point_def = 0;

	while (fgets(line, LINE_SIZE, f) != NULL)
	{
		token = strtok(line, " :");
	
		//Controllo parole chiave in file
		if (strncmp(token, "NAME", 4) == 0)
			continue;

		if (strncmp(token, "COMMENT", 7) == 0)
			continue;

		if (strncmp(token, "TYPE", 4) == 0)
			continue;

		if (strncmp(token, "DIMENSIONS", 10) == 0)
		{
			if(tsp_in->num_nodes<0)
				tsp_in->num_nodes=atoi(strtok(NULL, " :"));
	
			continue;
		}

		if (strncmp(token, "EDGE_WEIGHT_TYPE", 16) == 0)
			continue;

		if (strncmp(token, "NODE_COORD_SECTION", 18) == 0)
		{

			assert(tsp_in->num_nodes < 1); 			//Abbiamo già la dimensione del grafo?

			tsp_in->x_coords = (double*) calloc(tsp_in->num_nodes, sizeof(double));
			tsp_in->y_coords = (double*) calloc(tsp_in->num_nodes, sizeof(double));

			int count = 1;
			while (fgets(line, LINE_SIZE, f) != NULL && count<=(tsp_in->num_nodes))
			{
				int i = atoi(strtok(line, " "));
				
				assert(i<1 || i>(tsp_in->num_nodes)); //indice del nodo non valido

				tsp_in->x_coords[i - 1] = atof(strtok(NULL, " "));
				tsp_in->y_coords[i - 1] = atof(strtok(NULL, " "));

				count++;
			}

			assert(count < (tsp_in->num_nodes)); //meno nodi del numero dichiarato
			continue;
		}

		if (strncmp(token, "EOF", 3) == 0)
		{
			break;
		}
	}

	//VERBOSE
}