/**
	@file mtz_solver.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief MTZ solver.
*/

#include "mtz_solver.h"

void mtz_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in)
{
	mtz_build_model(env, lp, tsp_in);
	CPXmipopt(env, lp);

	int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
	tsp_in->sol = (double*)malloc(sizeof(double) * num_edges);

	assert(CPXgetmipx(env, lp, tsp_in->sol, 0, num_edges - 1) == 0);

	if (tsp_in->integerDist)
	{
		double cost;
		CPXgetbestobjval(env, lp, &cost);
		tsp_in->bestCostI = (int)(cost + CAST_PRECISION);
	}
	else
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
}

void mtz_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in)
{
#ifdef METAHEURISTIC
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, tsp_in->seed);

	printf(LINE);
	printf("%sMetaheuristic Procedure:%s \n", GREEN, WHITE);
	printf("%sSeed:%s %d\n\n", CYAN, WHITE,tsp_in->seed);

#endif

	//Properties of each edge
	char type = 'B';
	char** edge = calloc(sizeof(char*), 1);
	edge[0] = calloc(sizeof(char), NAME_SIZE);

	//i<j => j in [i+1, num_nodes-1]
	int i = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
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
			double ub = (i == j) ? 0.0 : 1.0;

			//Check if the column has been added correctly
			assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &type, edge) == 0);
			//Check if the index that we use to identify the edge
			//corresponds to the index that it should have
			assert(CPXgetnumcols(env, lp) - 1 == compact_xpos(tsp_in, i, j));
		}
	}

	type = 'I';
	i = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		sprintf(edge[0], "u%d", i + 1);

		double lb = 0.0;
		double ub = (i == 0) ? 0.0 : (double)(tsp_in->num_nodes - 2);

		double c = 0.0;

		//Check if the column has been added correctly
		assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &type, edge) == 0);
	}

	double const_term = 1.0;
	char type_constraint = 'E';
	char** constraint = calloc(sizeof(char*), 1);
	constraint[0] = calloc(sizeof(char), NAME_SIZE);

	int h = 0;
	for (; h < tsp_in->num_nodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		sprintf(constraint[0], "Degree-in constraint(%d)", h + 1);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)
		{
			assert(CPXchgcoef(env, lp, lastrow, compact_xpos(tsp_in, i, h), 1.0) == 0);
		}
	}

	for (h = 0; h < tsp_in->num_nodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		sprintf(constraint[0], "Degree-out constraint(%d)", h + 1);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{
			assert(CPXchgcoef(env, lp, lastrow, compact_xpos(tsp_in, h, j), 1.0) == 0);
		}
	}

	const_term = tsp_in->num_nodes - 1;
	type_constraint = 'L';
	int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
	int count = 0;
	i = 1;

#ifdef LAZY_CONSTRAINTS
	int nnz = 3;
	int indices[3];
	double values[3];
	int first = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			//u_i - u_j + n*x_ij <= n-1
			count++;
			sprintf(constraint[0], "Big-M constraint_lazy(%d)", count);
			indices[0] = compact_xpos(tsp_in, i, j);
			values[0] = tsp_in->num_nodes + 0.0;
			indices[1] = i + num_edges;
			values[1] = 1.0;
			indices[2] = j + num_edges;
			values[2] = -1.0;
			//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
			assert(CPXaddlazyconstraints(env, lp, 1, nnz, &const_term, &type_constraint, &first, indices, values, constraint) == 0); //one row for each node
		}
	}

	nnz = 2;
	int indices2[2];
	double values2[2];
	first = 0;
	const_term = 1;
	count = 0;
	i = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			//x_ij + x_ji <= 1
			count++;
			sprintf(constraint[0], "SEC 2nodes_constraint(%d)", count);
			indices[0] = compact_xpos(tsp_in, i, j);
			values[0] = 1.0;
			indices[1] = compact_xpos(tsp_in, j, i);
			values[1] = 1.0;
			//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
			assert(CPXaddlazyconstraints(env, lp, 1, nnz, &const_term, &type_constraint, &first, indices, values, constraint) == 0); //one row for each node
		}
	}

#else
	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			//u_i - u_j + n*x_ij <= n-1
			count++;
			int lastrow = CPXgetnumrows(env, lp);

			sprintf(constraint[0], "Big-M constraint(%d)", count);
			//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
			assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

			assert(CPXchgcoef(env, lp, lastrow, compact_xpos(tsp_in, i, j), tsp_in->num_nodes) == 0);

			if (i != j)
			{
				assert(CPXchgcoef(env, lp, lastrow, num_edges + i, 1.0) == 0);
				assert(CPXchgcoef(env, lp, lastrow, num_edges + j, -1.0) == 0);
			}
		}
	}
#endif

	if (tsp_in->verbose >= 100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(edge[0]);
	free(constraint[0]);
	free(edge);
	free(constraint);
}

void mtz_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp)
{
	int i = 1;
	for (; i < tsp_in->num_nodes; i++)
	{
		x[i] += 1.0;
		comp[i] = 1;
	}
	comp[0] = 1;

	double k = 1.0;
	int prev = 0;
	while (k < tsp_in->num_nodes - EPS)
	{
		for (i = 1; i < tsp_in->num_nodes; i++)
		{
			if (x[i] <= k + EPS && x[i] >= k - EPS)
			{
				succ[prev] = i;
				prev = i;
				k += 1.0;
			}
		}
	}
	succ[prev] = 0;
}
