#include "algorithms.h"

void default_alg(tsp_instance* tsp_in)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		tsp_in->sol[i] = i;
	}
	tsp_in->sol[tsp_in->num_nodes] = tsp_in->sol[0];
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
	
	//RICOMINCIARE DA QUI
	path_tsp(x, CPXgetnumcols(env, lp), tsp_in);
	
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

void path_tsp(double *x, int size, tsp_instance* tsp_in)
{
	int i = 0;
	int k = 0;

	double** matrix = calloc(sizeof(double), tsp_in->num_nodes);
	
	for(; i<tsp_in->num_nodes; i++)
		matrix[i]= calloc(sizeof(double), tsp_in->num_nodes);

	i = 0;

	for (; i < tsp_in->num_nodes && k<size; i++)
	{
		int j;
		for (j=i+1; j < tsp_in->num_nodes && k<size; j++)
		{
			matrix[i][j] = x[k];
			k++;
		}
	}

	int start = 0;
	int count = 0;
	int j;

	printf("%d \n",start+1);
	for (i = 0; i < tsp_in->num_nodes;)
	{
		int check = 1;
		
		for (j = i + 1; j < tsp_in->num_nodes && check; j++)
		{
			if (matrix[i][j] == 1.0)
			{
				matrix[i][j] = -1.0;
				printf("%d \n", j+1);
				count++;
				i = j;
				check = 0;
			}
		}

		if(check)
			j = i;
	
		for (i=0; i < j && check; i++)
		{
			if (matrix[i][j] == 1.0)
			{
				matrix[i][j] = -1.0;
				printf("%d \n", i+1);
				count++;
				check = 0;
			}
		}

		if (check)
		{
			start++;
			i = start;
			printf("start=%d \n", i+1);
		}
	}
}