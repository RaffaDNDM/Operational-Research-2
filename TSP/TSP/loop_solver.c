#include "loop_solver.h"

void loop_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp)
{
#ifdef METAHEURISTIC
	CPXsetintparam(env, CPX_PARAM_NODELIM, tsp_in->node_lim);
	CPXsetintparam(env, CPX_PARAM_POPULATELIM, tsp_in->sol_lim);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, tsp_in->eps_gap);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, tsp_in->seed);

	printf(LINE);
	printf("Metaheuristic Procedure: \n");
	printf("Node limit: %d    ", tsp_in->node_lim);
	printf("Pool solution: %d    ", tsp_in->sol_lim);
	printf("Cplex precision: %lf   ", tsp_in->eps_gap);
	printf("Seed: %d\n\n", tsp_in->seed);
#endif

	time_t start, start_iter, end_iter, end;
	start = clock();
	start_iter = start;
	cplex_build_model(tsp_in, env, lp);

	int i = 0;
	tsp_in->num_cols = CPXgetnumcols(env, lp);
	tsp_in->sol = calloc(sizeof(double), tsp_in->num_cols);
	//double* x = calloc(sizeof(double), tsp_in->num_cols);
	int n_comps = 3;

	CPXmipopt(env, lp);
	end_iter = clock();
	CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
	assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
	cplex_define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
	//cplex_plot(tsp_in, succ, comp, &n_comps);

	print_state(env, lp, n_comps, start_iter, end_iter);

	while (n_comps >= 2)
	{
		start_iter = clock();
		add_sec_constraint(env, lp, tsp_in, comp, n_comps);

		CPXmipopt(env, lp);
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
		end_iter = clock();
		assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
		cplex_define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
		//cplex_plot(tsp_in, succ, comp, &n_comps);
		print_state(env, lp, n_comps, start_iter, end_iter);

	}

	#ifdef METAHEURISTIC
		CPXsetintparam(env, CPX_PARAM_NODELIM, 9223372036800000000);
		CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 9223372036800000000);
		CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-04);

		start_iter = clock();
		CPXmipopt(env, lp);
		end_iter = clock();
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
		assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
		cplex_define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
		print_state(env, lp, n_comps, start_iter, end_iter);
		//cplex_plot(tsp_in, succ, comp, &n_comps);

		while (n_comps >= 2)
		{
			start_iter = clock();
			add_sec_constraint(env, lp, tsp_in, comp, n_comps);

			CPXmipopt(env, lp);
			end_iter = clock();
			CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
			assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
			cplex_define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
			print_state(env, lp, n_comps, start_iter, end_iter);
			//cplex_plot(tsp_in, succ, comp, &n_comps);
		}
	#endif

	//free(x);

	end = clock();

	if (tsp_in->integerDist)
		tsp_in->bestCostI = (int)(tsp_in->bestCostD + CAST_PRECISION);

	tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;
	
}


void add_sec_constraint(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* comp, int n_comps)
{
	static int iteration = 0;
	iteration++;

	int* const_terms = calloc(sizeof(int), n_comps);
	char type_constraint = 'L';
	char** constraint = calloc(sizeof(char*), 1);
	constraint[0] = calloc(sizeof(char), NAME_SIZE);

	int k = 0;
	for (; k < tsp_in->num_nodes; k++)
		const_terms[comp[k] - 1]++;

	for (k = 1; k <= n_comps; k++)
	{
		double const_term = const_terms[k - 1] - 1.0;
		sprintf(constraint[0], "SEC-constraint (%d, %d)", k, iteration);
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

		assert(CPXaddrows(env, lp, 0, 1, nnz, &const_term, &type_constraint, &first, indices, values, NULL, constraint) == 0); //one row for each node
	
		free(indices);
		free(values);
	}

	if (tsp_in->verbose >= 100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(constraint[0]);
	free(constraint);
	free(const_terms);
	

}

void print_state(CPXENVptr env, CPXLPptr lp, int ncomps, time_t start, time_t end)
{
	double actual_bound, best_bound;
	CPXgetbestobjval(env, lp, &best_bound);
	CPXgetobjval(env, lp, &actual_bound);
	printf("Best bound of all the remaining open nodes: %lf\n", best_bound);
	printf("Object value of the solution in the solution pool: %lf\n", actual_bound);
	printf("Number of components: %d\n", ncomps);
	printf("Time: %lf\n\n", (((double)(end - start) / (double)CLOCKS_PER_SEC)) * TIME_SCALE);
}