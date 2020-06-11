/**
	@file loop_solver.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Loop solver.
*/

#include "loop_solver.h"

void loop_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
{
#ifdef METAHEURISTIC
	CPXsetintparam(env, CPX_PARAM_NODELIM, tsp_in->node_lim);
	CPXsetintparam(env, CPX_PARAM_POPULATELIM, tsp_in->sol_lim);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, tsp_in->eps_gap);
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, tsp_in->seed);

	printf("%sMetaheuristic Procedure: %s\n", GREEN, WHITE);
	printf("%sNode limit:%s %d    ",CYAN, WHITE, tsp_in->node_lim);
	printf("%sPool solution:%s %d    ",CYAN, WHITE, tsp_in->sol_lim);
	printf("%sCplex precision:%s %lf   ",CYAN, WHITE, tsp_in->eps_gap);
	printf("%sSeed:%s %d\n\n",CYAN, WHITE, tsp_in->seed);
#endif

	time_t start, start_iter, end_iter, end;
	start = clock();
	start_iter = start;
	cplex_build_model(tsp_in, env, lp);

	int i = 0;
	tsp_in->num_cols = CPXgetnumcols(env, lp);
	tsp_in->sol = calloc(sizeof(double), tsp_in->num_cols);
	//double* x = calloc(sizeof(double), tsp_in->num_cols);

	CPXmipopt(env, lp);
	end_iter = clock();
	CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
	assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
	define_tour(tsp_in, tsp_in->sol, succ, comp, n_comps);
	//plot(tsp_in, succ, comp, &n_comps);

	print_state(env, lp, *n_comps, start_iter, end_iter);

	double remaining_time = tsp_in->deadline - ((double) (end_iter - start_iter) / CLOCKS_PER_SEC);

	while (*n_comps >= 2)
	{
		
		CPXsetdblparam(env, CPXPARAM_TimeLimit, remaining_time);
		start_iter = clock();
		add_sec_constraint(env, lp, tsp_in, comp, *n_comps);

		CPXmipopt(env, lp);
		end_iter = clock();
		remaining_time = remaining_time - ((double)(end_iter - start_iter) / CLOCKS_PER_SEC);
		if (remaining_time <= 0)
			break;
		
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
		assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
		//define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
		cplex_define_tour(tsp_in, tsp_in->sol, succ, comp, n_comps);
		
		print_state(env, lp, *n_comps, start_iter, end_iter);
		/*int j;
		for (j = 0; j < tsp_in->num_nodes; j++)
			printf("succ[%d] = %d    comp[%d] = %d \n", j, succ[j], j, comp[j]);*/

		//plot(tsp_in, succ, comp, &n_comps);

		
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
		define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
		print_state(env, lp, n_comps, start_iter, end_iter);
		//plot(tsp_in, succ, comp, &n_comps);

		remaining_time = tsp_in->deadline - ((double) (end_iter - start_iter) / CLOCKS_PER_SEC);

		while (*n_comps >= 2)
		{
			CPXsetdblparam(env, CPXPARAM_TimeLimit, remaining_time);
			start_iter = clock();
			add_sec_constraint(env, lp, tsp_in, comp, *n_comps);

			CPXmipopt(env, lp);
			end_iter = clock();
			CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
			assert(CPXgetmipx(env, lp, tsp_in->sol, 0, CPXgetnumcols(env, lp) - 1) == 0);
			define_tour(tsp_in, tsp_in->sol, succ, comp, n_comps);
			print_state(env, lp, *n_comps, start_iter, end_iter);

			remaining_time = remaining_time - ((double)(end_iter - start_iter) / CLOCKS_PER_SEC);
			if (remaining_time < 0)
				break;
			//plot(tsp_in, succ, comp, &n_comps);
		}
	#endif

	//free(x);

	//end = clock();

	if (tsp_in->integerDist)
	{
		tsp_in->bestCostI = (int)(tsp_in->bestCostD + CAST_PRECISION);
		tsp_in->bestCostD = 0.0;
	}

	//tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);

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

					indices[nnz] = xpos(tsp_in, i, j);
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
	printf("\n%sBest bound of all the remaining open nodes:%s %.2lf\n",BLUE,WHITE, best_bound);
	printf("%sObject value of the solution in the solution pool:%s %.2lf\n",BLUE, WHITE, actual_bound);
	printf("%sNumber of components:%s %d\n",BLUE, WHITE, ncomps);
	printf("%sTime:%s %.2lf\n\n",BLUE, WHITE, (((double)(end - start) / (double)CLOCKS_PER_SEC)) * TIME_SCALE);
}
