#include "cplex_solver.h"

void cplex_solver(tsp_instance* tsp_in)
{
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	cplex_build_model(tsp_in, env, lp);

	CPXmipopt(env, lp);

	if (tsp_in->integerDist)
	{
		double cost;
		CPXgetbestobjval(env, lp, &cost);
		tsp_in->bestCostI = (int)(cost + CAST_PRECISION);
	}
	else
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);

	//Values of variables in the solution
	double* x = calloc(sizeof(double), CPXgetnumcols(env, lp));
	assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);

	print_cost(tsp_in);

	if (tsp_in->plot)
	{
		int n_comps;
		int* succ = calloc(tsp_in->num_nodes, sizeof(int));
		int* comp = calloc(tsp_in->num_nodes, sizeof(int));

		define_tour(x, tsp_in, succ, comp, &n_comps);
		plot_cplex(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
	}

	free(x);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}

void cplex_build_model(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp)
{
	//Properties of each edge
	char bin = 'B';
	char** edge = calloc(sizeof(char*), 1);
	edge[0] = calloc(sizeof(char), NAME_SIZE);

	//i<j => j in [i+1, num_nodes-1]
	int i = 0;

	for (; i < (tsp_in->num_nodes - 1); i++)
	{
		int j;
		for (j = i + 1; j < tsp_in->num_nodes; j++)
		{
			sprintf(edge[0], "x%d,%d", i + 1, j + 1);

			double c; //cost of the edge

			if (tsp_in->integerDist)
			{
				int x;
				dist(i, j, tsp_in, &x);
				c = (double)x;
			}
			else
				dist(i, j, tsp_in, &c);

			double lb = 0.0;
			double ub = 1.0;

			//Check if the column has been added correctly
			assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &bin, edge) == 0);
			//Check if the index that we use to identify the edge 
			//corresponds to the index that it should have
			assert(CPXgetnumcols(env, lp) - 1 == xpos(tsp_in, i, j));
		}
	}


	double const_term = 2.0;
	char type_constraint = 'E';
	char** constraint = calloc(sizeof(char*), 1);
	constraint[0] = calloc(sizeof(char), NAME_SIZE);

	int h = 0;
	for (; h < tsp_in->num_nodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		sprintf(constraint[0], "Degree constraint(%d)", h + 1);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{
			if (j == h)
				continue;

			assert(CPXchgcoef(env, lp, lastrow, xpos(tsp_in, j, h), 1.0) == 0);
		}
	}

	if (tsp_in->verbose >= 100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(edge[0]);
	free(constraint[0]);
	free(edge);
	free(constraint);
}

int xpos(tsp_instance* tsp_in, int i, int j)
{
	assert(i != j); //error if ring edge

	if (i > j)
		return xpos(tsp_in, j, i);

	return (tsp_in->num_nodes * i + j) - ((i + 1) * (i + 2)) / 2;
}

void define_tour(double* x, tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
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

	#endif

	for (i = 0; i < tsp_in->num_nodes; succ[i++] = -1);

	*n_comps = 0;
	int begin;

	for (begin = 0; begin < tsp_in->num_nodes; begin++)
	{
		if (succ[begin] == -1)
		{
			(*n_comps)++;
			comp[begin] = n_comps;

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

void plot_cplex(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
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