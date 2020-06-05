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
		tsp_in->heu_sol = (int*)calloc((size_t)ncores, sizeof(int*));
		tsp_in->present_heu_sol = (int*)calloc((size_t)ncores, sizeof(int));; //non è ancora stata trovata una soluzione
		tsp_in->cost_heu_sol = (double*)calloc((size_t)ncores, sizeof(double));
		int k;
		for (k = 0; k < ncores; k++)
		{
			tsp_in->cost_heu_sol[k] = CPX_INFBOUND;
			tsp_in->heu_sol[k] = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		}

		//assert(CPXsetintparam(env, CPXPARAM_MIP_Strategy_HeuristicFreq, -1) == 0);

		CPXsetintparam(env, CPX_PARAM_LBHEUR, CPX_OFF);

		CPXsetlazyconstraintcallbackfunc(env, patching_callback, tsp_in);
		//CPXsetheuristiccallbackfunc(env, heuristic_callback, tsp_in);
		CPXsetincumbentcallbackfunc(env, set_incumbent, tsp_in);
		//CPXsetintparam(env, CPX_PARAM_HEURFREQ, 5);
		break;
	}
	}

	if (tsp_in->alg != 4)
		CPXsetintparam(env, CPX_PARAM_THREADS, ncores);

	int i = 0;
	int percentage[] = {90, 75, 50, 25, 10, 0};

	tsp_in->num_cols = CPXgetnumcols(env, lp);
	tsp_in->sol = (double*)calloc((size_t)tsp_in->num_cols, sizeof(double));

	int k = 2;
	int size_freedom = 10;
	int freedom[10];

	for (; k < 11; k++)
		freedom[k - 2] = k;

	freedom[k-2] = 20;

	int stop = 0;
	double cost = 0.0;
	tsp_in->bestCostD = CPX_INFBOUND;

	switch (tsp_in->heuristic)
	{
		case 1:
		{
			time_t start = 0;
			time_t end = 0;

			for (; i < 6; i++)
			{
				if (i == 0)
				{
					start = clock();
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_in->deadline);
				}
				else if (i == 1)
				{
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, DEFAULT_SOLLIM_VALUE);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, (tsp_in->deadline - ((double)(end - start) / (double)CLOCKS_PER_SEC)) / 5);
				}

				CPXmipopt(env, lp);
				end = clock();

				assert(CPXgetobjval(env, lp, &cost) == 0);

				if (tsp_in->bestCostD - cost > 1e-6)
				{
					tsp_in->bestCostD = cost;

					printf("-------- Improvement --------\n");

					assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);

				}

				cplex_change_coeff(tsp_in, env, lp, tsp_in->sol, percentage[i]);
			}
		}
		break;

		case 2:
		{
			tsp_in->bestCostD = CPX_INFBOUND;

			time_t start = 0;
			time_t end = 0;

			i = 0;
			for (;i< 10 && tsp_in->heuristic ; i++)
			{
				if (i == 0)
				{
					start = clock();
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, 1);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, tsp_in->deadline);
				}
				else if (i == 1)
				{
					CPXsetintparam(env, CPX_PARAM_INTSOLLIM, DEFAULT_SOLLIM_VALUE);
					CPXsetdblparam(env, CPXPARAM_TimeLimit, (tsp_in->deadline - ((double)(end - start) / (double)CLOCKS_PER_SEC)) / 9 );
				}

				CPXmipopt(env, lp);

				double gap;
				assert(CPXgetmiprelgap(env, lp, &gap) == 0);
				printf("\n gap = %lf \n", gap);

				end = clock();

				assert(CPXgetobjval(env, lp, &cost) == 0);

				if (tsp_in->bestCostD - cost > 1e-6)
				{
					tsp_in->bestCostD = cost;

					printf("-------- Improvement --------\n");

					assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);
				}
				local_branching(tsp_in, env, lp, tsp_in->sol, freedom[i]);
			}
		}
		break;
	}

	CPXmipopt(env, lp);

	/*
	if (tsp_in->alg == 4)//patching
	{
		int* succ = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		int* comp = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		int n_comps = 0;
		double* x = (double*)calloc((size_t)tsp_in->num_cols, sizeof(double));
		assert(CPXgetmipx(env, lp, x, 0, tsp_in->num_cols - 1) == 0);
		define_tour(tsp_in, x, succ, comp, &n_comps);
		
		while (n_comps != 1)
		{
			CPXmipopt(env, lp);
			define_tour(tsp_in, x, succ, comp, &n_comps);
		}

		free(succ);
		free(comp);
		free(x);
	}
	*/
	
	if (tsp_in->bestCostD - cost > 1e-6)
	{
		tsp_in->bestCostD = cost;

		printf("-------- Improvement --------\n");

		assert(CPXgetmipx(env, lp, tsp_in->sol, 0, tsp_in->num_cols - 1) == 0);
	}

	if (tsp_in->integerDist)
	{
		tsp_in->bestCostI = (int)(tsp_in->bestCostD + CAST_PRECISION);
		tsp_in->bestCostD = 0.0;
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

//controllare valore di ritorno, cosa restituire
static int CPXPUBLIC patching_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	//set quale thread deve lavorare

	double* x_star = (double*)malloc(tsp_in->num_cols * sizeof(double));

	if (CPXgetcallbacknodex(env, cbdata, wherefrom, x_star, 0, tsp_in->num_cols - 1))
	{
		free(x_star);
		return 1;
	}

	double objval;
	CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);

	patching(env, tsp_in, x_star, objval, cbdata, wherefrom);

	*useraction_p = CPX_CALLBACK_SET;

	free(x_star);

	return 0;
}

int patching(CPXENVptr env, tsp_instance* tsp_in, double* x_star, double objval, void* cbdata, int wherefrom)
{
	double cost = objval;

	int n_comps = 0;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	cplex_define_tour(tsp_in, x_star, succ, comp, &n_comps);

	//int h;
	//for (h = 0; h < tsp_in->num_nodes; h++)
	//printf("node %d -->> succ = %d -->> comp = %d\n", h, succ[h], comp[h]);

	//plot(tsp_in, succ, comp, &n_comps);

	//int num_comps_tot = n_comps;

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
			while( visited_nodes2[(i + 1) % num_nodes_vn2] != succ_node1)
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

		/*
		int k;
		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			if (comp[k] == n_comps - 1)
				printf("%d %d   n_comp: %d\n", k, succ[k], n_comps - 1);
		}
		printf("merged component %d\n", n_comps - 1);
		*/

		//printf("\ncomp2 = %d ------ \n", comp2);
		//printf("node1 = %d -->> succ[node1] = %d\n", node1, succ_node1);
		//printf("node2 = %d -->> succ[node2] = %d\n\n", node2, succ_node2);

		//for (i = 0; i < tsp_in->num_nodes; i++)
		//printf("node %d -->> succ = %d -->> comp = %d\n", i, succ[i], comp[i]);

		n_comps--;

		//for (i = 0; i < tsp_in->num_nodes; i++)
		//printf("i=%d    succ[i]=%d     comp[i]=%d\n", i, succ[i], comp[i]);

		//plot(tsp_in, succ, comp, &n_comps);

		free(visited_nodes1);
		free(visited_nodes2);
	}

	//plot(tsp_in, succ, comp, &n_comps);


	int thread;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &thread);
	printf(">>>>>>>>>>>>>>>>>>>>>>>>thread %d<<<<<<<<<<<<<<<<<<<<<<<<<<\n", thread);

	//if (cost < objval && cost < tsp_in->cost_heu_sol)
	if (cost < tsp_in->cost_heu_sol[thread])
	{
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

	if (tsp_in->present_heu_sol[thread] == 1)
	{


		printf("++++++++++++ Used heuristic solution ++++++++++++++   cost = %lf\n", tsp_in->cost_heu_sol[thread]);
		int i;
		for (i = 0; i < tsp_in->num_cols; i++)
			x[i] = 0.0;

		for (i = 0; i < tsp_in->num_nodes; i++)
			x[xpos(tsp_in, i, tsp_in->heu_sol[thread][i])] = 1.0;

		*objval_p = tsp_in->cost_heu_sol[thread];

		tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
		tsp_in->present_heu_sol[thread] = 0; //soluzione non più valida

		/*
		printf("thread %d\n", thread);
		printf("num cols %d\n", tsp_in->num_cols);

		for (i = 0; i < tsp_in->num_cols; i++)
		{
			if (x[i] == 1.0)
				printf("x[%d] = %lf\n", i, x[i]);
		}*/
		

		*useraction_p = CPX_CALLBACK_SET;
		*checkfeas_p = 1;
	}

	return 0;

}

static int CPXPUBLIC set_incumbent(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double objval, double* x, int* isfeas_p, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tsp_instance* tsp_in = (tsp_instance*) cbhandle;


	printf("objval = %lf\n", objval);

	int thread;
	CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, (void*)&thread);

	if (tsp_in->present_heu_sol[thread] == 1)
	{


		printf("++++++++++++ Used heuristic solution ++++++++++++++   cost = %lf\n", tsp_in->cost_heu_sol[thread]);
		int i;
		for (i = 0; i < tsp_in->num_cols; i++)
			x[i] = 0.0;

		for (i = 0; i < tsp_in->num_nodes; i++)
			x[xpos(tsp_in, i, tsp_in->heu_sol[thread][i])] = 1.0;

		//*objval_p = tsp_in->cost_heu_sol[thread];

		tsp_in->cost_heu_sol[thread] = CPX_INFBOUND;
		tsp_in->present_heu_sol[thread] = 0; //soluzione non più valida


		(*useraction_p) = CPX_CALLBACK_SET;
		(*isfeas_p) = 1;
	}


	return 0;
}
