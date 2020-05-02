#include "bc_solver.h"

void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general)
{
	cplex_build_model(tsp_in, env, lp);

	double epgap;
	CPXgetdblparam(env, CPX_PARAM_EPGAP, &epgap);
	printf("\n\n\nparam epgap: %lf\n\n\n", epgap);

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);

	if (general)
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, sec_general_callback, tsp_in);
	else
		CPXsetlazyconstraintcallbackfunc(env, sec_callback, tsp_in);

	int ncores = 1;
	CPXgetnumcores(env, &ncores);
	CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
	
	int i = 0;
	int percentage[] = {90, 75, 50, 25, 10, 0};

	tsp_in->num_cols = CPXgetnumcols(env, lp);
	double* x_best = malloc(sizeof(double) * tsp_in->num_cols);
	
	int k = 2;
	int freedom[10];
	for (; k < 11; k++)
		freedom[k - 2] = k;
	freedom[k-2] = 20;

	switch (tsp_in->heuristic)
	{
		case 1:
			for (; i < 6 ; i++)
			{
				if (i == 0)
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
				else if (i == 1)
				{
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 10000000);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, (tsp_in->deadline + 0.0) / 6.0);
				}

				CPXmipopt(env, lp);
				double cost;
				assert(CPXgetbestobjval(env, lp, &cost) == 0);
				assert(CPXgetmipx(env, lp, x_best, 0, tsp_in->num_cols - 1) == 0);
				cplex_change_coeff(tsp_in, env, lp, x_best, percentage[i]);
			}
			break;

		case 2:
		{
			double best_cost = CPX_INFBOUND;
			double cost;
			double cost2;
			double gap;
			k = 0;
			for (; i < 6 && tsp_in->heuristic; i++)
			{
				if (i == 0)
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
				else if (i == 1)
				{
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 10000000);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, (tsp_in->deadline + 0.0) / 6.0);
				}

				CPXmipopt(env, lp);
				assert(CPXgetmipx(env, lp, x_best, 0, tsp_in->num_cols - 1) == 0);
				assert(CPXgetobjval(env, lp, &cost) == 0);
				assert(CPXgetmiprelgap(env, lp, &gap) == 0);
				assert(CPXgetbestobjval(env, lp, &cost2) == 0);
				
				if (best_cost > cost)
				{
					k = 0;
					best_cost = cost;
					local_branching(tsp_in, env, lp, x_best, freedom[k]);
				}
				else
					local_branching(tsp_in, env, lp, x_best, freedom[++k]);

				printf("\n\n\n\ncost: %lf best cost2: %lf gap: %lf\n\n\n\n", cost,cost2, gap);
			}
		}
		break;
	}
	
	CPXmipopt(env, lp);

	free(x_best);
}

static int CPXPUBLIC sec_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	double* x_star = (double*)malloc(tsp_in->num_cols * sizeof(double));

	if (CPXgetcallbacknodex(env, cbdata, wherefrom, x_star, 0, tsp_in->num_cols - 1))
	{
		free(x_star);
		return 1;
	}

	int ncuts = sec_bc_constraint(env, tsp_in, x_star, cbdata, wherefrom);

	if (ncuts >= 1)
		*useraction_p = CPX_CALLBACK_SET;

	free(x_star);

	return 0;
}

static int CPXPUBLIC sec_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle)
{
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	double* x_star = (double*)malloc(tsp_in->num_cols * sizeof(double));

	if (CPXcallbackgetcandidatepoint(context, x_star, 0, tsp_in->num_cols - 1, NULL))
	{
		free(x_star);
		return 1;
	}

	int ncuts = sec_bc_constraint_general(context, tsp_in, x_star);
	free(x_star);

	return 0;
}

int sec_bc_constraint(CPXENVptr env, tsp_instance* tsp_in, double* x_star, void* cbdata, int wherefrom)
{
	int n_comps = 3;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));
	cplex_define_tour(tsp_in, x_star, succ, comp, &n_comps);

	if (n_comps == 1)
		return 0;

	int* const_terms = calloc(sizeof(int), n_comps);
	char type_constraint = 'L';

	int k = 0;
	for (; k < tsp_in->num_nodes; k++)
		const_terms[comp[k] - 1]++;

	for (k = 1; k <= n_comps; k++)
	{
		double const_term = const_terms[k - 1] - 1.0;
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		//assert(CPXnewrows(env, lp, 1, &(const_term), &type_constraint, NULL, constraint) == 0); //one row for each node

		int nnz = 0;
		int size = (tsp_in->num_nodes) / 4;
		int* indices = (int*)calloc(sizeof(int), size);
		double* values = (double*)calloc(sizeof(double), size);

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)
		{
			int j;
			for (j = i + 1; j < tsp_in->num_nodes; j++)
			{
				if (comp[i] == k && comp[j] == k)
				{
					if (nnz == size)
					{
						size *= 2;
						assert((indices = (int*)realloc(indices, size * sizeof(int))) != NULL);
						assert((values = (double*)realloc(values, size * sizeof(double))) != NULL);
					}

					indices[nnz] = cplex_xpos(tsp_in, i, j);
					values[nnz] = 1.0;
					nnz++;
				}
			}
		}

		if (nnz < size)
		{
			assert((indices = (int*)realloc(indices, nnz * sizeof(int))) != NULL);
			assert((values = (double*)realloc(values, nnz * sizeof(double))) != NULL);
		}

		assert(CPXcutcallbackadd(env, cbdata, wherefrom, nnz, const_term, type_constraint, indices, values, 0) == 0);

		free(indices);
		free(values);
	}

	free(const_terms);
	free(succ);
	free(comp);

	return n_comps;
}

int sec_bc_constraint_general(CPXCALLBACKCONTEXTptr context, tsp_instance* tsp_in, double* x_star)
{
	int n_comps = 3;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));
	cplex_define_tour(tsp_in, x_star, succ, comp, &n_comps);

	if (n_comps == 1)
		return 0;

	int* const_terms = calloc(sizeof(int), n_comps);
	char type_constraint = 'L';

	int k = 0;
	for (; k < tsp_in->num_nodes; k++)
		const_terms[comp[k] - 1]++;

	for (k = 1; k <= n_comps; k++)
	{
		double const_term = const_terms[k - 1] - 1.0;
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		//assert(CPXnewrows(env, lp, 1, &(const_term), &type_constraint, NULL, constraint) == 0); //one row for each node

		int nnz = 0;
		int size = (tsp_in->num_nodes) / 4;
		int* indices = (int*)calloc(sizeof(int), size);
		double* values = (double*)calloc(sizeof(double), size);
		int first = 0;

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)
		{
			int j;
			for (j = i + 1; j < tsp_in->num_nodes; j++)
			{
				if (comp[i] == k && comp[j] == k)
				{
					if (nnz == size)
					{
						size *= 2;
						assert((indices = (int*)realloc(indices, size * sizeof(int))) != NULL);
						assert((values = (double*)realloc(values, size * sizeof(double))) != NULL);
					}

					indices[nnz] = cplex_xpos(tsp_in, i, j);
					values[nnz] = 1.0;
					nnz++;
				}
			}
		}

		if (nnz < size)
		{
			assert((indices = (int*)realloc(indices, nnz * sizeof(int))) != NULL);
			assert((values = (double*)realloc(values, nnz * sizeof(double))) != NULL);
		}

		assert(CPXcallbackrejectcandidate(context, 1, nnz, &const_term, &type_constraint, &first, indices, values) == 0);

		free(indices);
		free(values);
	}

	free(const_terms);
	free(succ);
	free(comp);

	return n_comps;
}

void local_branching(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp, double* x_best, int freedom)
{
	static int call_num = 0;
	call_num++;

	printf("%s Iteration: %d freedom: %d\n%s", LINE, call_num, freedom, LINE);

	int i = 0;
	int nnz = 0;
	int first = 0;
	char type_constraint = 'G';
	double const_term = tsp_in->num_nodes - freedom + 1.0;
	int *indices = (int *) calloc(sizeof(int), tsp_in->num_nodes + 1);
	double *values = (double *) calloc(sizeof(double), tsp_in->num_nodes + 1);
	char** constraint = calloc(sizeof(char*), 1);
	constraint[0] = calloc(sizeof(char), NAME_SIZE);
	sprintf(constraint[0], "Local branching %d", call_num);

	if (call_num > 1)
		assert(CPXdelrows(env, lp, CPXgetnumrows(env, lp) - 1, CPXgetnumrows(env, lp) - 1)==0);

	for (; i < (tsp_in->num_nodes - 1) && nnz <= tsp_in->num_nodes; i++)
	{
		int j;
		for (j = i + 1; j < tsp_in->num_nodes && nnz <= tsp_in->num_nodes; j++)
		{
			int pos = cplex_xpos(tsp_in, i, j);

			assert(x_best[pos] <= EPS || x_best[pos] >= (1.0 - EPS));

			if (x_best[pos] > 0.5)
			{
				indices[nnz] = pos;
				values[nnz] = 1.0;
				nnz++;
			}
		}
	}
	
	assert(CPXaddrows(env, lp, 0, 1, nnz, &const_term, &type_constraint, &first, indices, values, NULL, constraint) == 0);

	CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(indices);
	free(values);
	free(constraint[0]);
	free(constraint);
}


void cplex_change_coeff(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp, double* x_best, int percentage)
{
	int i = 0;
	char which_bound = 'L';
	srand(time(NULL));
	double lb;

	for (; i < (tsp_in->num_nodes - 1); i++)
	{
		int j;
		for (j = i + 1; j < tsp_in->num_nodes; j++)
		{
			int pos = cplex_xpos(tsp_in, i, j);
			lb = 0.0;

			if (x_best[pos] > 0.5)
			{
				int choice = rand() % 100;
				lb = (choice < percentage) ? 1.0 : 0.0;
			}
			//printf("x (%d,%d) lb=%lf\n", i + 1, j + 1, lb);
			assert(CPXchgbds(env, lp, 1, &pos, &which_bound, &lb) == 0);
		}
	}
	CPXwriteprob(env, lp, LP_FILENAME, NULL);
}