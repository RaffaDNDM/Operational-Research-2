#include "utility.h"
#include <math.h>

void dist(int node1, int node2, tsp_instance* tsp_in, void* dist)
{
	double x_dist = tsp_in->x_coords[node1] - tsp_in->x_coords[node2];
	double y_dist = tsp_in->y_coords[node1] - tsp_in->y_coords[node2];

	if (tsp_in->integerDist)
	{
		int x_distI = (int)(x_dist + CAST_PRECISION);
		int y_distI = (int)(y_dist + CAST_PRECISION);

		int* distI = (int*)dist;
		*distI = (int)(sqrt((double)x_distI * (double)x_distI + (double)y_distI * (double)y_distI) + CAST_PRECISION);

		return;
	}

	double* distD = (double*)dist;
	*distD = sqrt(x_dist * x_dist + y_dist * y_dist);
}

void print_cost(tsp_instance* tsp_in)
{
	printf(LINE);
	printf("The solution cost is:\n");
	if (tsp_in->integerDist)
		printf("int cost %d\n", tsp_in->bestCostI);
	else
		printf("double cost %10.31f\n", tsp_in->bestCostD);
	printf(LINE);
}

void define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps)
{
	int i = 0;

#ifdef SOLUTION_CORRECTNESS

	int* degree = calloc(tsp_in->num_nodes, sizeof(int));
	for (; i < tsp_in->num_nodes; i++)
	{
		int j;
		for (j = i + 1; j < tsp_in->num_nodes; j++)
		{
			assert(x[xpos(tsp_in, i, j)] <= EPS || x[xpos(tsp_in, i, j)] >= (1.0 - EPS));

			if (x[xpos(tsp_in, i, j)] > 0.5)
			{
				degree[i]++;
				degree[j]++;
			}
		}
	}

	for (i = 0; i < tsp_in->num_nodes && degree[i] == 2; i++);
	assert(i == tsp_in->num_nodes);

	free(degree);

#endif

	for (i = 0; i < tsp_in->num_nodes; succ[i++] = -1)
		comp[i] = 0;

	*n_comps = 0;
	int begin;

	for (begin = 0; begin < tsp_in->num_nodes; begin++)
	{
		if (succ[begin] == -1)
		{
			(*n_comps)++;
			comp[begin] = *n_comps;

			int j = begin;

			for (i = 0; i < tsp_in->num_nodes; i++)
			{
				if (i == j)
					continue;

				if (x[xpos(tsp_in, j, i)] > 0.5 && comp[i] == 0)
				{
					succ[j] = i;
					comp[i] = *n_comps;
					j = i;
					i = -1;
				}
			}

			succ[j] = begin;
		}
	}
}

int xpos(tsp_instance* tsp_in, int i, int j)
{
	assert(i != j); //error if ring edge

	if (i > j)
		return xpos(tsp_in, j, i);

	return (tsp_in->num_nodes * i + j) - ((i + 1) * (i + 2)) / 2;
}

void plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
{
	FILE* f = fopen(CPLEX_DAT, "w");

	int count_comp = 1;
	for (; count_comp <= (*n_comps); count_comp++)
	{
		int i;
		for (i = 0; i < tsp_in->num_nodes && comp[i] != count_comp; i++);

		if (i == tsp_in->num_nodes)
			break;

		int begin = i; //first node of the count_comp-th component

		fprintf(f, "%f ", tsp_in->x_coords[begin]);
		fprintf(f, "%f ", tsp_in->y_coords[begin]);
		fprintf(f, "%d \n", begin + 1);

		int check = 1;

		int node = -1;
		do
		{
			node = succ[i];
			fprintf(f, "%f ", tsp_in->x_coords[node]);
			fprintf(f, "%f ", tsp_in->y_coords[node]);
			fprintf(f, "%d \n", node + 1);
			i = node;
		} while (node != begin);

		fprintf(f, "\n\n");
	}

	fclose(f);

	FILE* pipe = _popen(GNUPLOT_EXE, "w");
	f = fopen(GP_CPLEX_STYLE, "r");

	char line[LINE_SIZE];
	while (fgets(line, LINE_SIZE, f) != NULL)
	{
		if (strcmp(line, "LINE\n") == 0)
		{
			_pclose(pipe);
			printf("Type something to continue and create the image");
			gets(line, LINE_SIZE);
			pipe = _popen(GNUPLOT_EXE, "w");
			continue;
		}
		
		fprintf(pipe, "%s ", line);
	}

	_pclose(pipe);
	fclose(f);
}