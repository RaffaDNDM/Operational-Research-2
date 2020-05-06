#include "cplex_solver.h"
#include "bc_solver.h"
#include "gg_solver.h"
#include "loop_solver.h"
#include "mtz_solver.h"

void cplex_solver(tsp_instance* tsp_in)
{
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	tsp_in->bestCostI = -1;
	tsp_in->bestCostD = -1.0;
	tsp_in->sol = NULL;

	if (tsp_in->verbose > 60)
		CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);

	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	switch (tsp_in->alg)
	{
		case 1:
		{
			printf(LINE);
			printf("Loop solver\n\n");
			loop_solver(env, lp, tsp_in, succ, comp);
			break;
		}

		case 2:
		{
			printf(LINE);
			time_t start = clock();
			printf("Branch&Cut solver\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 0);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 3:
		{
			printf(LINE);
			time_t start = clock();
			printf("Branch&Cut solver with general callbacks\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 1);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 4:
		{
			printf(LINE);
			printf("MTZ solver\n\n");
			time_t start_time = clock();
			mtz_solver(env, lp, tsp_in);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 5:
		{
			printf(LINE);
			printf("GG solver\n\n");
			time_t start_time = clock();
			gg_solver(env, lp, tsp_in);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;

		}
	}
	
	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
			case 2:; case 3:
			{
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
				break;
			}

			case 4:
			{
				tsp_in->num_cols = CPXgetnumcols(env, lp);
				int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
				int start_pos = num_edges;
				int end_pos = tsp_in->num_cols - 1;

				double* x = (double*)malloc(sizeof(double) * tsp_in->num_nodes);
				assert(CPXgetmipx(env, lp, x, start_pos, end_pos) == 0);
			
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				mtz_define_tour(tsp_in, x, succ, comp);
				break;
			}

			case 5:
			{
				int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
				tsp_in->num_cols = CPXgetnumcols(env, lp);
				double* x = (double*)malloc(sizeof(double) * num_edges);

				int start_pos = num_edges;
				int end_pos = tsp_in->num_cols - 1;
				assert(CPXgetmipx(env, lp, x, start_pos, end_pos) == 0);

				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				gg_define_tour(tsp_in, x, succ, comp);
				break;
			}
		}

		plot(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
		free(tsp_in->sol);
	}

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

			double lb= 0.0;			
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

int cplex_xpos(tsp_instance* tsp_in, int i, int j)
{
	assert(i != j); //error if ring edge

	if (i > j)
		return cplex_xpos(tsp_in, j, i);

	return (tsp_in->num_nodes * i + j) - ((i + 1) * (i + 2)) / 2;
}

void cplex_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps)
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

void cplex_plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
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

int compact_xpos(tsp_instance* tsp_in, int i, int j)
{
	return i * (tsp_in->num_nodes) + j;
}
