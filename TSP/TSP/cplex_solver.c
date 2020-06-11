/**
	@file cplex_solver.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief CPLEX solvers.
*/

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

	if (tsp_in->deadline != DEADLINE_MAX)
		CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_in->deadline);

	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	int n_comps = 1;

	printf("%sCplex solver%s\n", RED, WHITE);
	printf("%s[Cplex]%s", BLUE, WHITE);

	switch (tsp_in->alg)
	{
		case 1:
		{
			time_t start = clock();
			printf("Loop solver\n\n");
			loop_solver(env, lp, tsp_in, succ, comp, &n_comps);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 2:
		{
			time_t start = clock();
			printf("Branch&Cut solver with general callbacks\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 0);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 3:
		{
			time_t start = clock();
			printf("Branch&Cut solver\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 1);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 4:
		{
			printf("Branch&Cut solver with patching\n\n");
			time_t start_time = clock();
			bc_solver(env, lp, tsp_in, succ, comp, 0);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 5:
		{
			printf("Branch&Cut solver with general callback and patching\n\n");
			time_t start_time = clock();
			bc_solver(env, lp, tsp_in, succ, comp, 0);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 6:
		{
			printf("MTZ solver\n\n");
			time_t start_time = clock();
			mtz_solver(env, lp, tsp_in);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;
		}

		case 7:
		{
			printf("GG solver\n\n");
			time_t start_time = clock();
			gg_solver(env, lp, tsp_in);
			time_t end_time = clock();
			tsp_in->execution_time = ((double)(end_time - start_time) / (double)CLOCKS_PER_SEC);
			break;

		}
	}

	print_cost(tsp_in);
	printf("%sExecution time:%s %.3lf seconds\n",GREEN, WHITE ,tsp_in->execution_time);
	printf("%s%s%s", RED, LINE, WHITE);

	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
		case 2:; case 3:; case 4:; case 5:
			{
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
				break;
			}

			case 6:
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

			case 7:
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

int compact_xpos(tsp_instance* tsp_in, int i, int j)
{
	return i * (tsp_in->num_nodes) + j;
}
