#include "algorithms.h"

void default_alg(tsp_instance* tsp_in)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		tsp_in->sol[i] = i;
	}
	tsp_in->sol[tsp_in->num_nodes] = tsp_in->sol[0];

	evaluate_sol(tsp_in);
	print_sol(tsp_in);
	
	if(tsp_in->plot)
		plot_solution(tsp_in);
}

void evaluate_sol(tsp_instance* tsp_in)
{
	int i = 1;
	int costI = 0;
	double costD = 0.0;

	tsp_in->bestCostD = 0;
	tsp_in->bestCostI = 0;

	if (tsp_in->integerDist)
	{
		int distI;

		for (; i < tsp_in->num_nodes + 1; i++)
		{
			dist(tsp_in->sol[i - 1], tsp_in->sol[i], tsp_in, &distI);
			costI += distI;
		}

		tsp_in->bestCostI = costI;
	}
	else
	{
		double distD;

		for (; i < tsp_in->num_nodes + 1; i++)
		{
			dist(tsp_in->sol[i - 1], tsp_in->sol[i], tsp_in, &distD);
			costD += distD;
		}

		tsp_in->bestCostD = costD;
	}
}

void cplex_solver(tsp_instance* tsp_in) 
{
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP");
	
	cplex_build_model(tsp_in, env, lp);

	CPXmipopt(env, lp);

	if (tsp_in->integerDist)
	{
		double cost;
		CPXgetbestobjval(env, lp, &cost);
		tsp_in->bestCostI = (int) (cost+CAST_PRECISION);
	}
	else
		CPXgetbestobjval(env, lp, &tsp_in->bestCostD);

	//Values of variables in the solution
	double* x = calloc(sizeof(double), CPXgetnumcols(env,lp));
	assert(CPXgetmipx(env, lp, x, 0, CPXgetnumcols(env, lp)-1)==0);
	
	print_sol(tsp_in);
	
	if (tsp_in->plot)
		plot_cplex(x, CPXgetnumcols(env, lp), tsp_in);
	
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}

void cplex_build_model(tsp_instance* tsp_in, CPXENVptr env, CPXLPptr lp)
{
	//Properties of each edge
	char bin = 'B';
	char** edge = calloc(sizeof(char*),1);
	edge[0] = calloc(sizeof(char),NAME_SIZE);

	//i<j => j in [i+1, num_nodes-1]
	int i = 0;

	for (; i < (tsp_in->num_nodes - 1); i++)
	{
		int j;
		for (j = i + 1; j < tsp_in->num_nodes; j++)
		{
			sprintf(edge[0], "x%d,%d", i+1, j+1);

			double c; //cost of the edge
			
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, j, tsp_in, &x);
				c = (double) x;
			}
			else
				dist(i, j, tsp_in, &c);

			double lb = 0.0;
			double ub = 1.0;

			//Check if the column has been added correctly
			assert(CPXnewcols(env, lp, 1, &c, &lb, &ub, &bin, edge) == 0);
			//Check if the index that we use to identify the edge 
			//corresponds to the index that it should have
			assert(CPXgetnumcols(env, lp) - 1 == xpos(tsp_in, i, j));
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
		sprintf(constraint[0], "Degree constraint(%d)",h+1);
		//CPXnewrows(env, lp, numero righe, vettore di termini noti, vettore di tipo di vincoli, NULL, cname)
		assert(CPXnewrows(env, lp, 1, &const_term, &type_constraint, NULL, constraint)==0); //one row for each node

		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{
			if (j == h)
				continue;

			assert(CPXchgcoef(env, lp, lastrow, xpos(tsp_in, j, h), 1.0)==0);
		}
	}

	if(tsp_in->verbose>=100)
		CPXwriteprob(env, lp, LP_FILENAME, NULL);

	free(edge[0]);
	free(constraint[0]);
	free(edge);
	free(constraint);
}

int xpos(tsp_instance* tsp_in, int i, int j)
{
	assert(i != j); //error if ring edge

	if (i > j)
		return xpos(tsp_in, j, i);

	return (tsp_in->num_nodes*i + j) - ((i+1)*(i+2))/2;
}

void plot_cplex(double *x, int size, tsp_instance* tsp_in)
{
	int i = 0;
	int k = 0;

	FILE* f = fopen(CPLEX_DAT, "w");

	for (; i < tsp_in->num_nodes && k<size; i++)
	{
		int j;
		for (j=i+1; j < tsp_in->num_nodes && k<size; j++)
		{
			if (x[k] == 1.0)
			{
				fprintf(f, "%f ", tsp_in->x_coords[i]);
				fprintf(f, "%f ", tsp_in->y_coords[i]);
				fprintf(f, "\n");
				fprintf(f, "%f ", tsp_in->x_coords[j]);
				fprintf(f, "%f ", tsp_in->y_coords[j]);
				fprintf(f, "\n\n");		
			}

			k++;
		}
	}

	fclose(f);

	FILE* pipe = _popen(GNUPLOT_EXE, "w");
	f = fopen(GP_CPLEX_STYLE, "r");

	char line[LINE_SIZE];
	while (fgets(line, LINE_SIZE, f) != NULL)
	{
		fprintf(pipe, "%s ", line);
	}

	fclose(f);
	_pclose(pipe);
}


void print_sol(tsp_instance* tsp_in)
{
	printf(LINE);
	printf("The solution cost is:\n");
	if (tsp_in->integerDist)
		printf("int cost %d\n", tsp_in->bestCostI);
	else
		printf("double cost %10.31f\n", tsp_in->bestCostD);
	printf(LINE);
}