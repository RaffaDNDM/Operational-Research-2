#include "cplex_solver.h"

void cplex_solver(tsp_instance* tsp_in)
{
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");

	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	switch (tsp_in->alg)
	{
		case 1:
		{
			printf(LINE);
			printf("Loop solver\n\n");
			loop_solver(env, lp, tsp_in, succ, comp);
			break;
		}

		case 2:
		{
			printf(LINE);
			time_t start = clock();
			printf("Branch&Cut solver\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 0);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;
			break;
		}

		case 3:
		{
			printf(LINE);
			time_t start = clock();
			printf("Branch&Cut solver with general callbacks\n\n");
			bc_solver(env, lp, tsp_in, succ, comp, 1);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;
			break;
		}

		case 4:
		{
			printf(LINE);
			printf("MTZ solver\n\n");
			time_t start = clock();
			mtz_build_model(env, lp, tsp_in);
			CPXmipopt(env, lp);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;
			break;
		}

		case 5:
		{
			printf(LINE);
			printf("GG solver\n\n");
			time_t start = clock();
			gg_build_model(env, lp, tsp_in);
			CPXmipopt(env, lp);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;
			break;

		}
	}

	if (tsp_in->integerDist)
	{
		double cost;
		CPXgetbestobjval(env, lp, &cost);
		tsp_in->bestCostI = (int)(cost + CAST_PRECISION);
	}
	else
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);

	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	//Values of variables in the solution
	double* x = calloc(sizeof(double), tsp_in->num_nodes);

	//switch su quante  variabili prendere per plottare

	switch (tsp_in->alg)
	{
		case 1:
		{
			x = (double*)realloc(x, sizeof(double) * CPXgetnumcols(env, lp));
			assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
			break;
		}

		case 2:; case 3:
		{
			x = (double*) realloc(x, sizeof(double) * CPXgetnumcols(env, lp));
			int s = CPXgetnumcols(env, lp);
			assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
			break;
		}

		case 4:
		{
			int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);

			int start = num_edges;
			int end = CPXgetnumcols(env, lp) - 1;

			assert(CPXgetmipx(env, lp, x, start, end ) == 0);
			break;
		}

		case 5:
		{
			int num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes);
			x = calloc(sizeof(double), num_edges);

			int start = num_edges;
			int end = CPXgetnumcols(env, lp) - 1;

			assert(CPXgetmipx(env, lp, x, start, end) == 0);
			break;
		}
	}

	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
			case 2:; case 3:
			{
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				cplex_define_tour(tsp_in, x, succ, comp, &n_comps);
				break;
			}

			case 4:
			{
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				mtz_define_tour(tsp_in, x, succ, comp);
				break;
			}

			case 5:
			{
				succ = calloc(tsp_in->num_nodes, sizeof(int));
				comp = calloc(tsp_in->num_nodes, sizeof(int));
				gg_define_tour(tsp_in, x, succ, comp);
				break;
			}
		}

		cplex_plot(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
	}

	free(x);
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

			double lb = 0.0;
			double ub = 1.0;

			//Check if the column has been added correctly
			assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &bin, edge) == 0);
			//Check if the index that we use to identify the edge
			//corresponds to the index that it should have
			assert(CPXgetnumcols(env, lp) - 1 == cplex_xpos(tsp_in, i, j));
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

			assert(CPXchgcoef(env, lp, lastrow, cplex_xpos(tsp_in, j, h), 1.0) == 0);
		}
	}

	if (tsp_in->verbose >= 100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(edge[0]);
	free(constraint[0]);
	free(edge);
	free(constraint);
}

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

static int CPXPUBLIC sec_callback(CPXCENVptr env, void* cbdata, int wherefrom, void* cbhandle, int* useraction_p)
{
	*useraction_p = CPX_CALLBACK_DEFAULT;
	tsp_instance* tsp_in = (tsp_instance*) cbhandle;

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
	for (;k < tsp_in->num_nodes; k++)
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
		for (;i < tsp_in->num_nodes; i++)
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
	for (;k < tsp_in->num_nodes; k++)
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
		for (;i < tsp_in->num_nodes; i++)
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

		assert(CPXcutcallbackadd(env, cbdata, wherefrom, nnz, const_term, type_constraint, indices, values, 0)==0);

		free(indices);
		free(values);
	}

	free(const_terms);
	free(succ);
	free(comp);

	return n_comps;
}

int cplex_xpos(tsp_instance* tsp_in, int i, int j)
{
	assert(i != j); //error if ring edge

	if (i > j)
		return cplex_xpos(tsp_in, j, i);

	return (tsp_in->num_nodes * i + j) - ((i + 1) * (i + 2)) / 2;
}

void cplex_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp, int* n_comps)
{
	int i = 0;

	#ifdef SOLUTION_CORRECTNESS

		int* degree = calloc(tsp_in->num_nodes, sizeof(int));
		for (; i < tsp_in->num_nodes; i++)
		{
			int j;
			for (j = i + 1; j < tsp_in->num_nodes; j++)
			{
				assert(x[cplex_xpos(tsp_in, i, j)] <= EPS || x[cplex_xpos(tsp_in, i, j)] >= (1.0 - EPS));

				if (x[cplex_xpos(tsp_in, i, j)] > 0.5)
				{
					degree[i]++;
					degree[j]++;
				}
			}
		}

		for (i = 0; i < tsp_in->num_nodes && degree[i] == 2; i++);
		assert(i == tsp_in->num_nodes);

	#endif

	for (i = 0; i < tsp_in->num_nodes; succ[i++] = -1)
		comp[i] = 0;

	*n_comps = 0;
	int begin;

	for (begin = 0; begin < tsp_in->num_nodes; begin++)
	{
		if (succ[begin] == -1)
		{
			(*n_comps)++;
			comp[begin] = *n_comps;

			int j = begin;

			for (i = 0; i < tsp_in->num_nodes; i++)
			{
				if (i == j)
					continue;

				if (x[cplex_xpos(tsp_in, j, i)] > 0.5 && comp[i] == 0)
				{
					succ[j] = i;
					comp[i] = *n_comps;
					j = i;
					i = -1;
				}
			}

			succ[j] = begin;
		}
	}
}

void cplex_plot(tsp_instance* tsp_in, int* succ, int* comp, int* n_comps)
{
	FILE* f = fopen(CPLEX_DAT, "w");

	int count_comp = 1;
	for (; count_comp <= (*n_comps); count_comp++)
	{
		int i;
		for (i = 0; i < tsp_in->num_nodes && comp[i] != count_comp; i++);

		if (i == tsp_in->num_nodes)
			break;

		int begin = i; //first node of the count_comp-th component

		fprintf(f, "%f ", tsp_in->x_coords[begin]);
		fprintf(f, "%f ", tsp_in->y_coords[begin]);
		fprintf(f, "%d \n", begin + 1);

		int check = 1;

		int node = -1;
		do
		{
			node = succ[i];
			fprintf(f, "%f ", tsp_in->x_coords[node]);
			fprintf(f, "%f ", tsp_in->y_coords[node]);
			fprintf(f, "%d \n", node + 1);
			i = node;
		} while (node != begin);

		fprintf(f, "\n\n");
	}

	fclose(f);

	FILE* pipe = _popen(GNUPLOT_EXE, "w");
	f = fopen(GP_CPLEX_STYLE, "r");

	char line[LINE_SIZE];
	while (fgets(line, LINE_SIZE, f) != NULL)
	{
		if (strcmp(line, "LINE\n") == 0)
		{
			_pclose(pipe);
			printf("Type something to continue and create the image");
			gets(line, LINE_SIZE);
			pipe = _popen(GNUPLOT_EXE, "w");
			continue;
		}

		fprintf(pipe, "%s ", line);
	}

	_pclose(pipe);
	fclose(f);
}

void mtz_build_model(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in)
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
			double ub = (i==j) ? 0.0 : 1.0;

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

	for (h=0; h < tsp_in->num_nodes; h++)
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

	const_term = tsp_in->num_nodes -1;
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
				values[0] = tsp_in->num_nodes+0.0;
				indices[1] = i+num_edges;
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

int compact_xpos(tsp_instance* tsp_in, int i, int j)
{
	return i * (tsp_in->num_nodes) + j;
}

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
			double ub = (i == j) ? 0.0 : ((j==0) ? 0.0: tsp_in->num_nodes-1);

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
		assert(CPXchgcoef(env, lp, lastrow, num_edges+compact_xpos(tsp_in, 0, j), 1.0) == 0);
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

void gg_define_tour(tsp_instance* tsp_in, double* x, int* succ, int* comp)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		comp[i] = 1;
	}

	double k = (double) (tsp_in->num_nodes-1);
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

void loop_solver(CPXENVptr env, CPXLPptr lp, tsp_instance* tsp_in, int* succ, int* comp)
{
	#ifdef METAHEURISTIC
		CPXsetintparam(env, CPX_PARAM_NODELIM, tsp_in->node_lim);
		CPXsetintparam(env, CPX_PARAM_POPULATELIM, tsp_in->sol_lim);
		CPXsetdblparam(env, CPX_PARAM_EPGAP, tsp_in->eps_gap);
		CPXsetintparam(env, CPX_PARAM_RANDOMSEED, tsp_in->seed);

		printf(LINE);
		printf("Metaheuristic Procedure: \n");
		printf("Node limit: %d    ", tsp_in->node_lim );
		printf("Pool solution: %d    ", tsp_in->sol_lim);
		printf("Cplex precision: %lf   ", tsp_in->eps_gap);
		printf("Seed: %d\n\n", tsp_in->seed);
	#endif

	time_t start, start_iter, end_iter, end;
	start = clock();
	start_iter = start;
	cplex_build_model(tsp_in, env, lp);

	int n_comps = 3;

	double* x = calloc(sizeof(double), CPXgetnumcols(env, lp));
	CPXmipopt(env, lp);
	end_iter = clock();
	CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
	assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
	cplex_define_tour(tsp_in, x, succ, comp, &n_comps);
	//cplex_plot(tsp_in, succ, comp, &n_comps);

	print_state(env, lp, n_comps, start_iter, end_iter);

	while (n_comps >= 2)
	{
		start_iter = clock();
		add_sec_constraint(env, lp, tsp_in, comp, n_comps);

		CPXmipopt(env, lp);
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
		end_iter = clock();
		assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
		cplex_define_tour(tsp_in, x, succ, comp, &n_comps);
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
		assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
		cplex_define_tour(tsp_in, x, succ, comp, &n_comps);
		print_state(env, lp, n_comps, start_iter, end_iter);
		//cplex_plot(tsp_in, succ, comp, &n_comps);

		while (n_comps >= 2)
		{
			start_iter = clock();
			add_sec_constraint(env, lp, tsp_in, comp, n_comps);

			CPXmipopt(env, lp);
			end_iter = clock();
			CPXgetbestobjval(env, lp, &tsp_in->bestCostD);
			assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp) - 1) == 0);
			cplex_define_tour(tsp_in, x, succ, comp, &n_comps);
			print_state(env, lp, n_comps, start_iter, end_iter);
			//cplex_plot(tsp_in, succ, comp, &n_comps);
		}
	#endif

	end = clock();

	tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC) * TIME_SCALE;

	free(x);
}

void add_sec_constraint(CPXENVptr env, CPXLPptr lp, tsp_instance *tsp_in, int *comp, int n_comps)
{
	static int iteration = 0;
	iteration++;

	int* const_terms = calloc(sizeof(int), n_comps);
	char type_constraint = 'L';
	char** constraint = calloc(sizeof(char*), 1);
	constraint[0] = calloc(sizeof(char), NAME_SIZE);

	int k = 0;
	for (;k < tsp_in->num_nodes; k++)
		const_terms[comp[k]-1]++;

	for (k=1; k <= n_comps; k++)
	{
		double const_term = const_terms[k-1]-1.0;
		sprintf(constraint[0], "SEC-constraint (%d, %d)", k, iteration);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		//assert(CPXnewrows(env, lp, 1, &(const_term), &type_constraint, NULL, constraint) == 0); //one row for each node

		int nnz = 0;
		int size = (tsp_in->num_nodes) / 4;
		int *indices = (int*) calloc(sizeof(int), size);
		double *values = (double*) calloc(sizeof(double), size);
		int first = 0;

		int i = 0;
		for (;i < tsp_in->num_nodes; i++)
		{
			int j;
			for (j = i + 1; j < tsp_in->num_nodes; j++)
			{
				if (comp[i] == k && comp[j] == k)
				{
					if (nnz == size)
					{
						size *= 2;
						assert((indices = (int*) realloc(indices, size * sizeof(int))) != NULL);
						assert((values = (double*) realloc(values, size * sizeof(double))) != NULL);
					}

					indices[nnz] = cplex_xpos(tsp_in, i, j);
					values[nnz] = 1.0;
					nnz++;
				}
			}
		}

		if (nnz < size)
		{
			assert((indices = (int*) realloc(indices, nnz * sizeof(int))) != NULL);
			assert((values = (double*) realloc(values, nnz * sizeof(double))) != NULL);
		}

		assert(CPXaddrows(env, lp, 0, 1, nnz, &const_term, &type_constraint, &first, indices, values, NULL, constraint) == 0); //one row for each node
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
	printf("Time: %lf\n\n",(((double)(end-start)/ (double) CLOCKS_PER_SEC)) * TIME_SCALE);
}
