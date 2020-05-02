#include "gg_solver.h"

void gg_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in)
{
#ifdef METAHEURISTIC
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, tsp_in->seed);

	printf(LINE);
	printf("Metaheuristic Procedure: \n");
	printf("Seed: %d\n\n", tsp_in->seed);

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



	//for yij => n*n variabili in pi�
	type = 'I';
	i = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{
			sprintf(edge[0], "y%d,%d", i + 1, j + 1);

			double c = 0.0; //cost of the edge

			double lb = 0.0;
			double ub = (i == j) ? 0.0 : ((j == 0) ? 0.0 : tsp_in->num_nodes - 1);

			//Check if the column has been added correctly
			assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &type, edge) == 0);
		}
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

	int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
	const_term = tsp_in->num_nodes - 1;
	type_constraint = 'E';

	int lastrow = CPXgetnumrows(env, lp);
	sprintf(constraint[0], "Flow-out source constraint");
	//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
	assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

	int j = 0;
	for (; j < tsp_in->num_nodes; j++)
	{
		assert(CPXchgcoef(env, lp, lastrow, num_edges + compact_xpos(tsp_in, 0, j), 1.0) == 0);
	}

	const_term = 1.0;
	type_constraint = 'E';

	h = 1;
	for (; h < tsp_in->num_nodes; h++)
	{
		int lastrow = CPXgetnumrows(env, lp);
		sprintf(constraint[0], "Balancing constraint(%d)", h + 1);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{
			assert(CPXchgcoef(env, lp, lastrow, num_edges + compact_xpos(tsp_in, j, h), 1.0) == 0);
			assert(CPXchgcoef(env, lp, lastrow, num_edges + compact_xpos(tsp_in, h, j), -1.0) == 0);
		}
	}

	const_term = 0.0;
	type_constraint = 'L';
	int coeff = -tsp_in->num_nodes + 1;

	i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		int j = 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			int lastrow = CPXgetnumrows(env, lp);
			sprintf(constraint[0], "Coupling constraint(%d,%d)", i + 1, j + 1);
			//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
			assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint) == 0); //one row for each node

			assert(CPXchgcoef(env, lp, lastrow, num_edges + compact_xpos(tsp_in, i, j), 1.0) == 0);
			assert(CPXchgcoef(env, lp, lastrow, compact_xpos(tsp_in, i, j), coeff) == 0);
		}
	}

	if (tsp_in->verbose >= 100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(edge[0]);
	free(constraint[0]);
	free(edge);
	free(constraint);
}

void gg_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		comp[i] = 1;
	}

	double k = (double)(tsp_in->num_nodes - 1);
	int prev = 0;

	i = 0;
	while (k > 0.0 + EPS)
	{
		int j = 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			int h = compact_xpos(tsp_in, i, j);

			if ((x[h] <= k + EPS) && (x[h] >= k - EPS))
			{
				succ[i] = j;
				i = j;
				break;
			}
		}
		k -= 1.0;
	}
	succ[i] = 0;
}
