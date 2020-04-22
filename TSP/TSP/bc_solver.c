#include "bc_solver.h"


void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general)
{
	cplex_build_model(tsp_in, env, lp);
	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	tsp_in->num_cols = CPXgetnumcols(env, lp);

	if (general)
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, sec_general_callback, tsp_in);
	else
		CPXsetlazyconstraintcallbackfunc(env, sec_callback, tsp_in);

	int ncores = 1;
	CPXgetnumcores(env, &ncores);
	CPXsetintparam(env, CPX_PARAM_THREADS, ncores);

	CPXmipopt(env, lp);
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