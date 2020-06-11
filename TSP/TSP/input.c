/**
	@file input.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Input management of TSPlib file and command line arguments.
*/

#include "input.h"

void parse_cmd(char** argv, int argc, tsp_instance* tsp_in)
{
	//at leat one param + program name
	if (argc == 1)
	{
		help();
	}
	//char** commands = malloc(sizeof(char*)*NUM_COMMANDS);

	tsp_in->num_nodes = -1;
	tsp_in->deadline = DEADLINE_MAX;
	tsp_in->alg = -1;
	tsp_in->integerDist = 0;
	tsp_in->plot = 1;
	tsp_in->verbose = VERBOSE;
	tsp_in->size = -1;
	tsp_in->heuristic = 0;

	int def_deadline = 0;

	strcpy(tsp_in->input, "NULL");
	strcpy(tsp_in->dir, "NULL");

	int i = 0;
	for (; i < NUM_ALGS; i++)
		tsp_in->which_alg[i] = 0;

	i = 1;
	for (; i < argc; i++)
	{
		if (strncmp(argv[i], "-alg", 4) == 0)
		{
			select_alg(tsp_in, argv[++i], 0);
			continue;
		}

		if (strncmp(argv[i], "-size", 5) == 0)
		{
			double sizeF = atof(argv[++i]);
			int size = (int)sizeF;

			//the value inserted by the user must be an integer (alg!=0 && algF==alg)
			//but also the value must be greater than zero
			assert(size > 0 && sizeF == size);

			tsp_in->size = size;

			continue;
		}

		if (strncmp(argv[i], "-int", 4) == 0 || strncmp(argv[i], "-i", 2) == 0)
		{
			tsp_in->integerDist = 1;
			continue;
		}

		if (strncmp(argv[i], "-v", 2) == 0 || strncmp(argv[i], "-verbose", 8) == 0)
		{
			double verbF = atof(argv[++i]);
			int verb = (int)verbF;

			//the value inserted by the user must be an integer (n_heu!=0 && heuF==n_heu)
			//but also the value must be greater than zero
			if (!(verb > 0 && verbF == verb))
			{
				printf("Not valid verbose value\n");
				return;
			}

			tsp_in->verbose = verb;
			continue;
		}

		if (strncmp(argv[i], "-noplot", 7) == 0 || strncmp(argv[i], "-np", 3) == 0)
		{
			tsp_in->plot = 0;
			continue;
		}

		if ((strncmp(argv[i], "-dead",5) == 0 || strncmp(argv[i], "-deadline",9) == 0) && !def_deadline)
		{
			double deadline =(double) atof(argv[++i]);

			//the value inserted by the user must be a floating point number (deadline!=0 )
			//but also the value must be greater than zero
			assert(deadline > 0 && deadline <= DEADLINE_MAX);

			def_deadline = 1;
			tsp_in->deadline = deadline;
			continue;
		}

		if ((strncmp(argv[i], "-f", 2) == 0 || strncmp(argv[i], "-file", 5) == 0 || strncmp(argv[i], "-input", 6) == 0))
		{
			if (strncmp(tsp_in->input, "NULL", 4) == 0)
			{
				strcpy(tsp_in->input, argv[++i]);
				continue;
			}

		}

		if (strncmp(argv[i], "-dir", 4) == 0)
		{
			if (strncmp(tsp_in->dir, "NULL", 4) == 0)
			{
				strcpy(tsp_in->dir, argv[++i]);
				continue;
			}

		}

		if (strncmp(argv[i], "-heu", 4) == 0)
		{
			double heuF = atof(argv[++i]);
			int n_heu = (int)heuF;

			//the value inserted by the user must be an integer (n_heu!=0 && heuF==n_heu)
			//but also the value must be greater than zero
			if (!(n_heu > 0 && heuF == n_heu))
			{
				printf("Undefined heuristic algorithm\n");
				return;
			}
			assert(n_heu> 0 && n_heu <= 2 );

			tsp_in->heuristic = n_heu;
			continue;
		}

		if ((strncmp(argv[i], "-help", 5) == 0 || strncmp(argv[i], "-h", 2) == 0))
		{
			//print set of commands and exit from the program
			help();
			exit(0);
		}
	}

	if (tsp_in->verbose > 30)
	{
		printf("%s%s%s",RED, LINE, WHITE);
		printf("%sList of parameters specified on command line:%s \n", RED,WHITE);
		int i;
		for (i=1; i<argc; i++)
		{
			int i_check = strncmp(argv[i], "-i", 2) == 0 || strncmp(argv[i], "-int", 4) == 0;
			int np_check = strncmp(argv[i], "-noplot", 7) == 0 || strncmp(argv[i], "-np", 3) == 0;

			if (i_check || np_check )
			{
				printf("%s%s%s\n",YELLOW, argv[i], WHITE);
			}
			else
			{
				printf("%s%s%s : %s\n",YELLOW, argv[i], WHITE, argv[i+1]);
				i++;
			}
		}
	}

	if (tsp_in->heuristic != 0)
	{
		if (tsp_in->alg != 2 && tsp_in->alg != 3)
		{
			printf("Not valid algorithm for heuristic option\n");
			exit(1);
		}
		else if (tsp_in->which_alg[0] == 1 || tsp_in->which_alg[3] == 1 || tsp_in->which_alg[4] == 1 ||
			tsp_in->which_alg[5] == 1 || tsp_in->which_alg[6] == 1 || tsp_in->which_alg[7] == 1 ||
			tsp_in->which_alg[8] == 1 || tsp_in->which_alg[9] == 1 || tsp_in->which_alg[10] == 1 ||
			tsp_in->which_alg[11] == 1)
		{
			printf("The heuristic option will be ignored for algorithms different from 2 and 3\n");
		}
	}
}

void help()
{
	printf(LINE);
	printf("                                       Help\n");
	printf(LINE);
	printf("Insert the algs you want to test\n");
	printf("-alg alg_num                where alg_num=string of numbers of algoritms (E.g. \"1235\")\n");
	printf("                                          (or \"all\" to use all the algorithms)\n\n");
	printf("Possible algorithms\n");
	printf("1) %s \n", ALG1);
	printf("2) %s \n", ALG2);
	printf("3) %s \n", ALG3);
	printf("4) %s \n", ALG4);
	printf("5) %s \n", ALG5);
	printf("6) %s \n", ALG6);
	printf("7) %s \n", ALG7);
	printf("8) %s \n", ALG8);
	printf("9) %s \n", ALG9);
	printf("10) %s \n", ALG10);
	printf("11) %s \n", ALG11);
	printf(STAR_LINE);
	printf("Insert the max time of the execution\n");
	printf("-d dead_time\n");
	printf("-dead dead_time             where dead_time = max execution time in seconds (float)\n");
	printf("-deadline dead_time\n");
	printf(STAR_LINE);
	printf("Insert the directory where there are input files\n");
	printf("-dir path_name              where path_name = existing path\n");
	printf(STAR_LINE);
	printf("Insert the file in input\n");
	printf("-f file_name.tsp\n");
	printf("-file file_name.tsp         where file_name = name of tsp file (input instance)\n");
	printf("-input file_name.tsp\n");
	printf(STAR_LINE);
	printf("Use integer costs\n");
	printf("-i\n");
	printf("-int\n");
	printf(STAR_LINE);
	printf("Set no plotting of solution \n");
	printf("-noplot\n");
	printf("-np\n");
	printf(STAR_LINE);
	printf("Insert the max number of nodes in instances in specified dir (with -dir)\n");
	printf("-size max_size              where max_size = max number of nodes for an instance\n");
	printf(STAR_LINE);
	printf("Set verbose information during the execution\n");
	printf("-v\n");
	printf("-verbose\n");
	printf(LINE);
	exit(0);
}


void parse_file(tsp_instance* tsp_in)
{
	FILE* f = fopen(tsp_in->input, "r");

	assert(f != NULL);

	tsp_in->num_nodes = -1;

	char line[LINE_SIZE];
	char* token;
	int point_def = 0;

	while (fgets(line, LINE_SIZE, f) != NULL)
	{
		token = strtok(line, " :");

		//Check key words in tsp file
		if (strncmp(token, "NAME", 4) == 0)
			continue;

		if (strncmp(token, "COMMENT", 7) == 0)
			continue;

		if (strncmp(token, "TYPE", 4) == 0)
			continue;

		if (strncmp(token, "DIMENSION", 9) == 0)
		{
			if (tsp_in->num_nodes < 0)
				tsp_in->num_nodes = atoi(strtok(NULL, " :"));

			continue;
		}

		if (strncmp(token, "EDGE_WEIGHT_TYPE", 16) == 0)
			continue;

		if (strncmp(token, "NODE_COORD_SECTION", 18) == 0)
		{
			//The number of nodes must be defined before this section
			assert(tsp_in->num_nodes > 0);

			tsp_in->x_coords = (double*) calloc(tsp_in->num_nodes, sizeof(double));
			tsp_in->y_coords = (double*) calloc(tsp_in->num_nodes, sizeof(double));

			int count = 1;
			while (fgets(line, LINE_SIZE, f) != NULL && count<=(tsp_in->num_nodes))
			{
				int i = atoi(strtok(line, " "));

				assert(i>0 && i<=(tsp_in->num_nodes)); //Index in {1,num_nodes}

				tsp_in->x_coords[i - 1] = atof(strtok(NULL, " "));
				tsp_in->y_coords[i - 1] = atof(strtok(NULL, " "));

				count++;
			}

			assert(count >= (tsp_in->num_nodes)); //Few nodes declarations
			continue;
		}

		if (strncmp(token, "EOF", 3) == 0)
		{
			break;
		}
	}

	fclose(f);

	printf("%sName of the input instance :%s %s\n",GREEN, WHITE, tsp_in->input);
	printf("%sNumber of input nodes :%s %d\n",GREEN, WHITE, tsp_in->num_nodes);

	if (tsp_in->verbose > 80)
	{
		if (tsp_in->deadline < DEADLINE_MAX)
			printf("%sDeadline time :%s %f\n",GREEN, WHITE, tsp_in->deadline);

		printf("\n%sInput nodes coordinates:%s\n", GREEN, WHITE);

		for (int i = 0; i < tsp_in->num_nodes; i++)
		{
			printf("node %3d : x = %10.2f  y=%10.2f\n", i + 1, tsp_in->x_coords[i], tsp_in->y_coords[i]);
		}
	}

	printf("%s%s%s", RED, LINE, WHITE);
}

void select_alg(tsp_instance* tsp_in, char* alg_string, int in_main)
{
	if (strncmp(alg_string, "all", 3) == 0)
	{
		tsp_in->alg = 0;

		int j = 0;
		for (; j < NUM_ALGS; j++)
		{
			tsp_in->which_alg[j] = 1;
		}
	}
	else
	{
		double algF = atof(alg_string);
		int alg = (int)algF;

		//the value inserted by the user must be an integer (alg!=0 && algF==alg)
		//but also the value must be greater than zero
		if (in_main)
			if (!(alg > 0 && algF == alg))
			{
				printf("Undefined algorithm\n");
				return;
			}
		else
			assert(alg > 0 && algF == alg);

		if (alg <= NUM_ALGS)
			tsp_in->alg = alg;
		else
		{
			tsp_in->alg = 0;
			int length = strlen(alg_string);

			int j = 0;
			for (; j < length; j++)
			{
				int index = (alg_string[j] - '0') - 1;

				if (index >= 0 && index < NUM_ALGS)
					tsp_in->which_alg[index] = 1;
				else
				{
					printf("Undefined algorithm\n");

					if (in_main)
					{
						tsp_in->alg = -1;
						return;
					}
					else
						exit(1);
				}
			}
		}
	}
}

void dealloc_inst(tsp_instance* tsp_in)
{
	free(tsp_in->x_coords);
	free(tsp_in->y_coords);
	//free(tsp_in->sol);
	tsp_in->x_coords = NULL;
	tsp_in->y_coords = NULL;
	tsp_in->sol = NULL;
}
