/**
	@file bc_solver.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief CPLEX solutions with Branch & Cut.
*/

#include "bc_solver.h"

void bc_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp, int general)
{
	cplex_build_model(tsp_in, env, lp);

	CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF);
	tsp_in->num_cols = CPXgetnumcols(env, lp);

	int ncores = 1;
	CPXgetnumcores(env, &ncores);

	switch (tsp_in->alg)
	{
	case 2:
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, sec_general_callback, tsp_in);
		break;

	case 3:
		CPXsetlazyconstraintcallbackfunc(env, sec_callback, tsp_in);
		break;

	case 4:
	{
		tsp_in->heu_sol = (int**)calloc((size_t)ncores, sizeof(int*));
		tsp_in->present_heu_sol = (int*)calloc((size_t)ncores, sizeof(int));; //non è ancora stata trovata una soluzione
		tsp_in->cost_heu_sol = (double*)calloc((size_t)ncores, sizeof(double));
		int k;
		for (k = 0; k < ncores; k++)
		{
			tsp_in->cost_heu_sol[k] = CPX_INFBOUND;
			tsp_in->heu_sol[k] = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		}

		CPXsetlazyconstraintcallbackfunc(env, sec_callback, tsp_in);
		CPXsetheuristiccallbackfunc(env, heuristic_callback, tsp_in);

		break;
	}
	case 5:
	{
		tsp_in->heu_sol = (int**)calloc((size_t)ncores, sizeof(int*));
		tsp_in->present_heu_sol = (int*)calloc((size_t)ncores, sizeof(int));; //non è ancora stata trovata una soluzione
		tsp_in->cost_heu_sol = (double*)calloc((size_t)ncores, sizeof(double));
		int k;
		for (k = 0; k < ncores; k++)
		{
			tsp_in->cost_heu_sol[k] = CPX_INFBOUND;
			tsp_in->heu_sol[k] = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		}

		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE | CPX_CALLBACKCONTEXT_LOCAL_PROGRESS, sec_general_callback, tsp_in);
		
		break;
	}
	}

	
	CPXsetintparam(env, CPX_PARAM_THREADS, ncores);

	int i = 0;
	int percentage[] = { 90, 75, 50, 25, 10, 0 };

	tsp_in->num_cols = CPXgetnumcols(env, lp);
	tsp_in->sol = (double*)calloc((size_t)tsp_in->num_cols, sizeof(double));

	int k = 2;
	int size_freedom = 10;
	int freedom[10];

	for (; k < 11; k++)
		freedom[k - 2] = k;

	freedom[k - 2] = 20;

	int stop = 0;
	double cost = 0.0;
	tsp_in->bestCostD = CPX_INFBOUND;

	switch (tsp_in->heuristic)
	{
	case 1:
	{
		time_t start = 0;
		time_t end = 0;

		double remaining_time = tsp_in->deadline;

		double step;

		int count = 0;

		for (; remaining_time >0 ; i++)
		{
			start = clock();

			if (i == 0)
			{
				CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
				CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_in->deadline);
			}
			else 
			{
				CPXsetintparam(env, CPX_PARAM_INTSOLLIM, DEFAULT_SOLLIM_VALUE);

				if (remaining_time < step)
					step = remaining_time;
				CPXsetdblparam(env, CPXPARAM_TimeLimit, step);
				//CPXsetdblparam(env, CPXPARAM_TimeLimit, (tsp_in->deadline - ((double)(end - start) / (double)CLOCKS_PER_SEC)) / 5);
			}

			CPXmipopt(env, lp);
			end = clock();

			if (i == 0)
				step = (tsp_in->deadline - ((double)(end - start) / (double)CLOCKS_PER_SEC)) / 5;


			remaining_time = remaining_time - ((double)(end - start) / (double)CLOCKS_PER_SEC);

			assert(CPXgetobjval(env, lp, &cost) == 0);

			if (tsp_in->bestCostD - cost > 1e-6)
			{
				tsp_in->bestCostD = cost;

				printf("%s-------- Improvement --------%s\n", GREEN, WHITE);

				assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);

				count = 0;

			}
			else
			{
				count++;
				if (count == 6)
					break;
			}

			cplex_change_coeff(tsp_in, env, lp, tsp_in->sol, percentage[i % 5]);
		}
	}
	break;

	case 2:
	{
		tsp_in->bestCostD = CPX_INFBOUND;

		time_t start = 0;
		time_t end = 0;

		double remaining_time = tsp_in->deadline;

		double step;
		int count = 0;

		i = 0;
		for (; remaining_time > 0; i++)
		{
			start = clock();

			if (i == 0)
			{
				
				CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
				CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_in->deadline);
			}
			else if (i != 0)
			{
				CPXsetintparam(env, CPX_PARAM_INTSOLLIM, DEFAULT_SOLLIM_VALUE);

				if (remaining_time < step)
					step = remaining_time;
				CPXsetdblparam(env, CPXPARAM_TimeLimit, step);
				//CPXsetdblparam(env, CPXPARAM_TimeLimit, remaining_time/9.0);
			}

			CPXmipopt(env, lp);

			end = clock();
			if(i==0)
				step =( tsp_in->deadline - ((double)(end - start) / (double)CLOCKS_PER_SEC)) / 9;

			remaining_time = remaining_time - ((double)(end - start) / (double)CLOCKS_PER_SEC);
				
			assert(CPXgetobjval(env, lp, &cost) == 0);
		
			if (tsp_in->bestCostD - cost > 1e-6)
			{
				tsp_in->bestCostD = cost;

				printf("%s-------- Improvement --------%s\n", GREEN, WHITE);

				assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);

				count = 0;
			}
			else
			{
				count++;
				if (count == 10)
					break;
			}
			local_branching(tsp_in, env, lp, tsp_in->sol, freedom[i % 10]);
		}
	}
	break;
	}

	CPXmipopt(env, lp);

	assert(CPXgetobjval(env, lp, &cost) == 0);

	if (tsp_in->bestCostD - cost > 1e-6)
	{
		tsp_in->bestCostD = cost;

		printf("%s-------- Improvement --------%s\n", GREEN, WHITE);

		assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);
	}

	if (tsp_in->integerDist)
	{
		tsp_in->bestCostI = (int)(tsp_in->bestCostD + CAST_PRECISION);
		tsp_in->bestCostD = 0.0;
	}

	if (tsp_in->alg == 4)
	{
		for (i=0; i < ncores; i++)
			free(tsp_in->heu_sol[i]);

		free(tsp_in->heu_sol);
		free(tsp_in->cost_heu_sol);
		free(tsp_in->present_heu_sol);
	}

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

	if (tsp_in->alg == 4)
	{
		if (CPXgetcallbacknodex(env, cbdata, wherefrom, x_star, 0, tsp_in->num_cols - 1))
		{
			free(x_star);
			return 1;
		}

		double objval;
		CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);

		int thread;
		CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thread);

		patching(tsp_in, x_star, objval, thread);
	}

	free(x_star);

	return 0;
}

static int CPXPUBLIC sec_general_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle)
{
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	int thread;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &thread);

	double* x_star = (double*)malloc(tsp_in->num_cols * sizeof(double));

	if (contextid == CPX_CALLBACKCONTEXT_CANDIDATE)
	{
		if (CPXcallbackgetcandidatepoint(context, x_star, 0, tsp_in->num_cols - 1, NULL))
		{
			free(x_star);
			return 1;
		}

		int ncuts = sec_bc_constraint_general(context, tsp_in, x_star);

		if (tsp_in->alg == 5)
		{
			double objval;

			if (CPXcallbackgetcandidatepoint(context, x_star, 0, tsp_in->num_cols - 1, &objval))
			{
				free(x_star);
				return 1;
			}

			patching(tsp_in, x_star, objval, thread);
		}
	}
	else if (tsp_in->alg == 5 && contextid == CPX_CALLBACKCONTEXT_LOCAL_PROGRESS)
	{
		if (thread == 0)
		{
			int thread_min = -1;
			int min_cost = CPX_INFBOUND;

			int ncores = 1;
			CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADS, &thread);

			int i;
			for (i = 0; i < ncores; i++)
			{
				//if (tsp_in->present_heu_sol[i] == 1)
					//printf("thread %d ----- cost %lf\n", i, tsp_in->cost_heu_sol[i]);

				if (tsp_in->present_heu_sol[i] == 1 && tsp_in->cost_heu_sol[i] < min_cost)
				{
					thread_min = i;
					min_cost = tsp_in->cost_heu_sol[i];
				}
			}

			if (thread_min != -1)
			{
				printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n",
					GREEN, YELLOW, tsp_in->cost_heu_sol[thread_min], RED, thread_min, GREEN, WHITE);

				int* indices = (int*)calloc((size_t)tsp_in->num_cols, sizeof(int));
				double* values = (double*)calloc((size_t)tsp_in->num_cols, sizeof(double));

				int i;
				for (i = 0; i < tsp_in->num_cols; i++)
					indices[i] = i;

				for (i = 0; i < tsp_in->num_nodes; i++)
					values[xpos(tsp_in, i, tsp_in->heu_sol[thread_min][i])] = 1.0;

				assert(CPXcallbackpostheursoln(context, tsp_in->num_cols, indices, values, tsp_in->cost_heu_sol[thread_min], CPXCALLBACKSOLUTION_NOCHECK) == 0);

				tsp_in->cost_heu_sol[thread_min] = CPX_INFBOUND;
				tsp_in->present_heu_sol[thread_min] = 0;

				free(indices);
				free(values);
			}
		}
		else
		{
			if (tsp_in->present_heu_sol[thread] == 1)
			{
				printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n",
					GREEN, YELLOW, tsp_in->cost_heu_sol[thread], RED, thread, GREEN, WHITE);

				int* indices = (int*)calloc((size_t)tsp_in->num_cols, sizeof(int));
				double* values = (double*)calloc((size_t)tsp_in->num_cols, sizeof(double));

				int i;
				for (i = 0; i < tsp_in->num_cols; i++)
					indices[i] = i;

				for (i = 0; i < tsp_in->num_nodes; i++)
					values[xpos(tsp_in, i, tsp_in->heu_sol[thread][i])] = 1.0;

				assert(CPXcallbackpostheursoln(context, tsp_in->num_cols, indices, values, tsp_in->cost_heu_sol[thread], CPXCALLBACKSOLUTION_NOCHECK) == 0);

				tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
				tsp_in->present_heu_sol[thread] = 0;

				free(indices);
				free(values);
			}
		}
		
		
	}

	free(x_star);

	return 0;
}

int sec_bc_constraint(CPXENVptr env, tsp_instance* tsp_in, double* x_star, void* cbdata, int wherefrom)
{
	int n_comps = 3;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));
	define_tour(tsp_in, x_star, succ, comp, &n_comps);

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
	define_tour(tsp_in, x_star, succ, comp, &n_comps);

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

	printf("%s%s%s", GREEN, LINE, WHITE);
	printf("%sIteration:%s %d %sfreedom:%s %d\n",CYAN, WHITE, call_num, CYAN, WHITE, freedom);
	printf("%s%s%s", GREEN, LINE, WHITE);

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
			int pos = xpos(tsp_in, i, j);

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
			int pos = xpos(tsp_in, i, j);
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

int patching(tsp_instance* tsp_in, double* x_star, double objval, int thread)
{
	double cost = objval;

	int n_comps = 0;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	cplex_define_tour(tsp_in, x_star, succ, comp, &n_comps);

	while (n_comps != 1)
	{

		double delta_min = CPX_INFBOUND;
		int node1 = -1;
		int succ_node1 = -1;
		int node2 = -1;
		int succ_node2 = -1;
		//int comp1 = -1;
		int comp2 = -1;

		int i;
		for (i = 0; i < tsp_in->num_nodes; i++)
		{
			if (comp[i] != 1)
				continue;

			int j;
			for (j = i+1; j < tsp_in->num_nodes; j++)
			{
				double cost_i_h; //[comp1[i], comp1[i+1]]
				double cost_j_k; //[comp2[j], comp2[j+1]]
				double cost_i_j; //[comp1[i], comp2[j]]
				double cost_h_k; //[comp1[i+1], comp2[j+1]]
				double cost_h_j; //[comp1[i+1], comp2[j]]
				double cost_k_i; //[comp2[j+1], comp1[i]]

				if (comp[j] != 1) //ho selezionato la prima componente , ora la espando unendola alla più vicina
				{

					if (tsp_in->integerDist)
					{
						int x1, x2, x3, x4, x5, x6;
						dist(i, succ[i], tsp_in, &x1);
						cost_i_h = (double)x1;
						dist(j, succ[j], tsp_in, &x2);
						cost_j_k = (double)x2;
						dist(i, j, tsp_in, &x3);
						cost_i_j = (double)x3;
						dist(succ[i], succ[j], tsp_in, &x4);
						cost_h_k = (double)x4;
						dist(succ[i], j, tsp_in, &x5);
						cost_h_j = (double)x5;
						dist(succ[j], i, tsp_in, &x6);
						cost_k_i = (double)x6;

					}
					else
					{
						dist(i, succ[i], tsp_in, &cost_i_h);
						dist(j, succ[j], tsp_in, &cost_j_k);
						dist(i, j, tsp_in, &cost_i_j);
						dist(succ[i], succ[j], tsp_in, &cost_h_k);
						dist(succ[i], j, tsp_in, &cost_h_j);
						dist(succ[j], i, tsp_in, &cost_k_i);
					}

					double delta1 = cost_i_j + cost_h_k - cost_i_h - cost_j_k;
					double delta2 = cost_h_j + cost_k_i - cost_i_h - cost_j_k;

					if (delta1 < delta2 && delta1 < delta_min)
					{
						delta_min = delta1;
						node1 = i;
						succ_node1 = j;
						node2 = succ[j];
						succ_node2 = succ[i];
						//comp1 = comp[i];
						comp2 = comp[j];
					}
					else if (delta2 < delta1 && delta2 < delta_min)
					{
						delta_min = delta2;
						node1 = i;
						succ_node1 = succ[j];
						node2 = j;
						succ_node2 = succ[i];
						//comp1 = comp[i];
						comp2 = comp[j];
					}
				}
			}
		}

		int* visited_nodes1 =(int*) calloc(tsp_in->num_nodes, sizeof(int));
		int* visited_nodes2 =(int*) calloc(tsp_in->num_nodes, sizeof(int));
		int num_nodes_vn1 = 0;
		int num_nodes_vn2 = 0;
		
		for (i = 0; comp[i] != 1; i++);

		int begin = i;
		visited_nodes1[num_nodes_vn1] = i;
		num_nodes_vn1++;

		while (succ[i] != begin)
		{
			//printf("%d %d \n", succ[i], comp[succ[i]]);

			visited_nodes1[num_nodes_vn1] = succ[i];
			num_nodes_vn1++;
			i = succ[i];
		}

		visited_nodes1 = (int*)realloc(visited_nodes1, num_nodes_vn1 * sizeof(int));

		for (i=0; comp[i] != comp2; i++);

		begin = i;
		visited_nodes2[num_nodes_vn2] = i;
		num_nodes_vn2++;

		while (succ[i] != begin)
		{
			visited_nodes2[num_nodes_vn2] = succ[i];
			num_nodes_vn2++;
			i = succ[i];
		}

		visited_nodes2 = (int*)realloc(visited_nodes2, num_nodes_vn2 * sizeof(int));

		cost = cost + delta_min;

		int change_order = 0;

		succ[node1] = succ_node1;
		succ[node2] = succ_node2;

		if (succ[succ_node1] == node2)
			change_order = 1;

		//printf("node1: %d ->>> succ_node1: %d");

		int count = 0;

		if (change_order)
		{
			for (i = 0; visited_nodes2[i] != node2; i++);

			//for (count=0; count < num_nodes_vn2 - 1; count++)
			//while( visited_nodes2[(i + 1) % num_nodes_vn2] != succ_node1)
			while (visited_nodes2[(i ) % num_nodes_vn2] != succ_node1)
			{
				succ[visited_nodes2[(i + 1) % num_nodes_vn2]] = visited_nodes2[i];
				i = (i + 1) % num_nodes_vn2;
			}
		}

		for (count = 0; count < tsp_in->num_nodes; count++)
		{
			if (comp[count] == comp2)
				comp[count] = 1;
			else if (comp[count] > comp2)
				comp[count]--;
		}

		n_comps--;

		//for (i = 0; i < tsp_in->num_nodes; i++)
			//printf("i=%d    succ[i]=%d     comp[i]=%d\n", i, succ[i], comp[i]);

		//printf("n_comps ---> %d\n", n_comps);
		//plot(tsp_in, succ, comp, &n_comps);

		free(visited_nodes1);
		free(visited_nodes2);
	}

	//plot(tsp_in, succ, comp, &n_comps);

	//if (cost < objval && cost < tsp_in->cost_heu_sol)
	if (cost < tsp_in->cost_heu_sol[thread])
	{
		printf("%scost heuristic solution%s %.2lf %s [thread %d]%s\n",BLUE, YELLOW, cost, RED,thread, WHITE);

		int i;
		for (i = 0; i < tsp_in->num_nodes; i++)
			tsp_in->heu_sol[thread][i] = 0;

		for (i = 0; i < tsp_in->num_nodes; i++)
			tsp_in->heu_sol[thread][i] = succ[i];

		tsp_in->cost_heu_sol[thread] = cost;

		tsp_in->present_heu_sol[thread] = 1;
	}

	free(succ);
	free(comp);
}

static int CPXPUBLIC heuristic_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;


	int thread;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM,(void *) &thread);
	/*
	if (tsp_in->present_heu_sol[thread] == 1)
	{
		printf("++++++++++++ Found heuristic solution [thread %d] ++++++++++++++   cost = %lf\n", thread, tsp_in->cost_heu_sol[thread]);
		
		int i;
		for (i = 0; i < tsp_in->num_cols; i++)
			x[i] = 0.0;

		for (i = 0; i < tsp_in->num_nodes; i++)
			x[xpos(tsp_in, i, tsp_in->heu_sol[thread][i])] = 1.0;

		*objval_p = tsp_in->cost_heu_sol[thread];

		tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
		tsp_in->present_heu_sol[thread] = 0; //soluzione non più valida

		*useraction_p = CPX_CALLBACK_SET;
		*checkfeas_p = 1;
	}*/
	
	if (thread == 0)
	{
		int thread_min = -1;
		int min_cost = CPX_INFBOUND;

		int ncores = 1;
		CPXgetnumcores(env, &ncores);
		
		int i;
		for (i = 0; i < ncores; i++)
		{
			//if (tsp_in->present_heu_sol[i] == 1)
				//printf("thread %d ----- cost %lf\n", i, tsp_in->cost_heu_sol[i]);

			if (tsp_in->present_heu_sol[i] == 1 && tsp_in->cost_heu_sol[i] < min_cost)
			{
				thread_min = i;
				min_cost = tsp_in->cost_heu_sol[i];
			}
		}
		//printf("\n");

		if (thread_min != -1)
		{
			printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n", GREEN, YELLOW, tsp_in->cost_heu_sol[thread_min],RED, thread, GREEN, WHITE);

			for (i = 0; i < tsp_in->num_cols; i++)
				x[i] = 0.0;

			for (i = 0; i < tsp_in->num_nodes; i++)
				x[xpos(tsp_in, i, tsp_in->heu_sol[thread_min][i])] = 1.0;

			*objval_p = tsp_in->cost_heu_sol[thread_min];

			tsp_in->cost_heu_sol[thread_min] = CPX_INFBOUND;
			tsp_in->present_heu_sol[thread_min] = 0; //soluzione non più valida

			*useraction_p = CPX_CALLBACK_SET;
			*checkfeas_p = 1;
		}
	}
	else
	{
		
		if (tsp_in->present_heu_sol[thread] == 1)
		{
			printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n", GREEN, YELLOW, tsp_in->cost_heu_sol[thread], RED, thread, GREEN, WHITE);

			int i;
			for (i = 0; i < tsp_in->num_cols; i++)
				x[i] = 0.0;

			for (i = 0; i < tsp_in->num_nodes; i++)
				x[xpos(tsp_in, i, tsp_in->heu_sol[thread][i])] = 1.0;

			*objval_p = tsp_in->cost_heu_sol[thread];

			tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
			tsp_in->present_heu_sol[thread] = 0; //soluzione non più valida

			*useraction_p = CPX_CALLBACK_SET;
			*checkfeas_p = 1;
		}
	}
	return 0;
}

static int CPXPUBLIC general_heuristic_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* cbhandle)
{
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	int thread;
	CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &thread);
	
	if (thread == 0)
	{
		int thread_min = -1;
		int min_cost = CPX_INFBOUND;

		int ncores = 1;
		CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADS, &thread);

		int i;
		for (i = 0; i < ncores; i++)
		{
			//if (tsp_in->present_heu_sol[i] == 1)
				//printf("thread %d ----- cost %lf\n", i, tsp_in->cost_heu_sol[i]);

			if (tsp_in->present_heu_sol[i] == 1 && tsp_in->cost_heu_sol[i] < min_cost)
			{
				thread_min = i;
				min_cost = tsp_in->cost_heu_sol[i];
			}
		}
		//printf("\n");

		if (thread_min != -1)
		{
			printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n", 
				GREEN, YELLOW, tsp_in->cost_heu_sol[thread_min], RED, thread, GREEN, WHITE);

			int* indices = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
			double* values = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));

			for (i = 0; i < tsp_in->num_nodes; i++)
			{
				indices[i] = xpos(tsp_in, i, tsp_in->heu_sol[thread_min][i]);
				values[i] = 1.0;
			}
		
			assert(CPXcallbackpostheursoln(context, tsp_in->num_nodes, indices, values, tsp_in->cost_heu_sol[thread_min], CPXCALLBACKSOLUTION_NOCHECK) == 0);

			tsp_in->cost_heu_sol[thread_min] = CPX_INFBOUND;
			tsp_in->present_heu_sol[thread_min] = 0; //soluzione non più valida
		}
	}
	else
	{
		if (tsp_in->present_heu_sol[thread] == 1)
		{
			printf("%s++++++++++++ Found heuristic solution with cost:%s %.2lf %s [thread %d]%s ++++++++++++++  %s\n", 
				GREEN, YELLOW, tsp_in->cost_heu_sol[thread], RED, thread, GREEN, WHITE);

			int i;
			int* indices = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
			double* values = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));

			for (i = 0; i < tsp_in->num_nodes; i++)
			{
				indices[i] = xpos(tsp_in, i, tsp_in->heu_sol[thread][i]);
				values[i] = 1.0;
			}

			assert(CPXcallbackpostheursoln(context, tsp_in->num_nodes, indices, values, tsp_in->cost_heu_sol[thread], CPXCALLBACKSOLUTION_NOCHECK) == 0);

			tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
			tsp_in->present_heu_sol[thread] = 0; //soluzione non più valida

		}
	}

	return 0;
}
