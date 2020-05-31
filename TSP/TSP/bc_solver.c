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

	switch (tsp_in->alg)
	{
	case 2:
		CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_CANDIDATE, sec_general_callback, tsp_in);
		break;

	case 3:
		CPXsetlazyconstraintcallbackfunc(env, sec_callback, tsp_in);
		break;

	case 8:
	{
		tsp_in->heu_sol = (int*)calloc((size_t)tsp_in->num_cols, sizeof(int));
		tsp_in->heu_sol[0] = -1; //non è ancora stata trovata una soluzione
		tsp_in->cost_heu_sol = CPX_INFBOUND;
		CPXsetlazyconstraintcallbackfunc(env, patching_callback, tsp_in);
		CPXsetheuristiccallbackfunc(env, heuristic_callback, tsp_in);
		break;
	}
	}

	if (tsp_in->alg != 8)
	{
		int ncores = 1;
		CPXgetnumcores(env, &ncores);
		CPXsetintparam(env, CPX_PARAM_THREADS, ncores);
	}
	else 
		CPXsetintparam(env, CPX_PARAM_THREADS, 1);

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
	double cost;
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
	assert(CPXgetobjval(env, lp, &cost) == 0);

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
}

int patching(CPXENVptr env, tsp_instance* tsp_in, double* x_star, double objval, void* cbdata, int wherefrom)
{
	double cost = objval;

	int n_comps = 0;
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	cplex_define_tour(tsp_in, x_star, succ, comp, &n_comps);
	//plot(tsp_in, succ, comp, &n_comps);

	while (n_comps != 1)
	{
		int* visited_nodes1 = calloc(tsp_in->num_nodes, sizeof(int));
		int* visited_nodes2 = calloc(tsp_in->num_nodes, sizeof(int));
		int num_nodes_vn1 = 0;
		int num_nodes_vn2 = 0;

		int k;
		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			if (comp[k] == n_comps)
				printf("%d %d   n_comp: %d\n", k, succ[k], n_comps);
		}
		printf("finish component %d\n", n_comps);

		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			if (comp[k] == n_comps - 1)
				printf("%d %d   n_comp: %d\n", k, succ[k], n_comps - 1);
		}
		printf("finish component %d\n", n_comps - 1);

		plot(tsp_in, succ, comp, &n_comps);


		int i = 0;

		for (; comp[i] != n_comps; i++);

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

		i = 0;
		for (; comp[i] != n_comps - 1; i++);

		begin = i;
		visited_nodes2[num_nodes_vn2++] = i;

		while (succ[i] != begin)
		{
			visited_nodes2[num_nodes_vn2++] = succ[i];
			i = succ[i];
		}

		visited_nodes2 = (int*)realloc(visited_nodes2, num_nodes_vn2 * sizeof(int));

		double delta_min = CPX_INFBOUND;
		int node1 = -1;
		int succ_node1 = -1;
		int node2 = -1;
		int succ_node2 = -1;

		for (i = 0; i < num_nodes_vn1 - 1; i++)
		{
			int j = 0;
			for (; j < num_nodes_vn2 - 1; j++)
			{
				double cost_i_h; //[comp1[i], comp1[i+1]]
				double cost_j_k; //[comp2[j], comp2[j+1]]
				double cost_i_j; //[comp1[i], comp2[j]]
				double cost_h_k; //[comp1[i+1], comp2[j+1]]
				double cost_h_j; //[comp1[i+1], comp2[j]]
				double cost_k_i; //[comp2[j+1], comp1[i]]

				if (tsp_in->integerDist)
				{
					int x1, x2, x3, x4, x5, x6;
					dist(visited_nodes1[i], visited_nodes1[i + 1], tsp_in, &x1);
					cost_i_h = (double)x1;
					dist(visited_nodes2[j], visited_nodes2[j + 1], tsp_in, &x2);
					cost_j_k = (double)x2;
					dist(visited_nodes1[i], visited_nodes2[j], tsp_in, &x3);
					cost_i_j = (double)x3;
					dist(visited_nodes1[i + 1], visited_nodes2[j + 1], tsp_in, &x4);
					cost_h_k = (double)x4;
					dist(visited_nodes1[i + 1], visited_nodes2[j], tsp_in, &x5);
					cost_h_j = (double)x5;
					dist(visited_nodes2[j + 1], visited_nodes1[i], tsp_in, &x6);
					cost_k_i = (double)x6;

				}
				else
				{
					dist(visited_nodes1[i], visited_nodes1[i + 1], tsp_in, &cost_i_h);
					dist(visited_nodes2[j], visited_nodes2[j + 1], tsp_in, &cost_j_k);
					dist(visited_nodes1[i], visited_nodes2[j], tsp_in, &cost_i_j);
					dist(visited_nodes1[i + 1], visited_nodes2[j + 1], tsp_in, &cost_h_k);
					dist(visited_nodes1[i + 1], visited_nodes2[j], tsp_in, &cost_h_j);
					dist(visited_nodes2[j + 1], visited_nodes1[i], tsp_in, &cost_k_i);
				}

				double delta1 = cost_i_h + cost_j_k - cost_i_j - cost_h_k;
				double delta2 = cost_i_h + cost_j_k - cost_h_j - cost_k_i;

				if (delta1 < delta2 && delta1 < delta_min)
				{
					delta_min = delta1;
					node1 = visited_nodes1[i];
					succ_node1 = visited_nodes2[j];
					node2 = visited_nodes1[i + 1];
					succ_node2 = visited_nodes2[j + 1];
				}
				else if (delta2 < delta1 && delta2 < delta_min)
				{
					delta_min = delta2;
					node1 = visited_nodes1[i + 1];
					succ_node1 = visited_nodes2[j];
					node2 = visited_nodes2[j + 1];
					succ_node2 = visited_nodes1[i];
				}

			}
		}

		cost = cost + delta_min;

		int change_order = 0;

		if (succ[node1] != succ_node2 || succ[node2] != succ_node1)
		{
			//change_order = 1;
			succ[succ_node1] = node1;
			succ[succ_node2] = node2;

		}
		else
		{
			succ[node1] = succ_node1;
			succ[node2] = succ_node2;

		}

		int count = 0;
		for (; count < tsp_in->num_nodes; count++)
		{
			if (succ[count] == succ_node1 || succ[count] == node1)
			{
				change_order = 1;
				break;
			}
		}


		//if (succ[succ_node1] == node2 || succ[succ_node1] == succ_node2 || succ[succ_node2] == node1 || succ[succ_node2] == succ_node1)
		if (change_order)
		{

			if (comp[succ_node1] == n_comps)
				for (i = 0; visited_nodes1[i] != succ_node1; i++);
			else
				for (i = 0; visited_nodes1[i] != succ_node2; i++);

			//devo invertire l'ordine, scansione modulare circolare, n-1 iterazioni

			for (count=0; count < num_nodes_vn1 - 1; count++)
			{
				succ[visited_nodes1[(i + 1) % num_nodes_vn1]] = visited_nodes1[i];
				//succ[visited_nodes1[i]] = visited_nodes1[(i - 1) % num_nodes_vn1];
				comp[visited_nodes1[(i + 1) % num_nodes_vn1]] = n_comps - 1;
				i = (i + 1) % num_nodes_vn1;
			}

			comp[visited_nodes1[i]] = n_comps - 1;
		}
		else
		{
			/*
			if (comp[node2] == n_comps)
				i = node2;
			else
				i = succ_node2;

			for (; count < num_nodes_vn1; count++)
			{
				i = succ[i];
				comp[i]=n_comps - 1;
			}*/

			for (count = 0; count < tsp_in->num_nodes; count++)
			{
				if (comp[count] == n_comps)
					comp[count]--;
			}
		}

		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			if (comp[k] == n_comps - 1)
				printf("%d %d   n_comp: %d\n", k, succ[k], n_comps - 1);
		}
		printf("merged component %d\n", n_comps - 1);

		n_comps--;

		plot(tsp_in, succ, comp, &n_comps);

		printf("\n\ndone\n");

		free(visited_nodes1);
		free(visited_nodes2);
	}

	if (cost < objval && cost < tsp_in->cost_heu_sol)
	{
		int i;
		for (i = 0; i < tsp_in->num_cols; i++)
			tsp_in->heu_sol[i] = 0.0;

		for (i = 0; i < tsp_in->num_nodes; i++)
			tsp_in->heu_sol[xpos(tsp_in, i, succ[i])] = 1.1;

		tsp_in->cost_heu_sol = cost;
	}

	free(succ);
	free(comp);
}

static int CPXPUBLIC heuristic_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, double* objval_p, double* x, int* checkfeas_p, int* useraction_p)
{
	tsp_instance* tsp_in = (tsp_instance*)cbhandle;

	if (tsp_in->heu_sol[0] != -1)
	{
		int i;
		for (i = 0; i < tsp_in->num_cols; i++)
			x[i] = tsp_in->heu_sol[i];

		*objval_p = tsp_in->cost_heu_sol;

		tsp_in->cost_heu_sol = CPX_INFBOUND;
		tsp_in->heu_sol[0] = -1.0; //soluzione non più valida
	}
}