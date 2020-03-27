#include "default_alg.h"

void default_alg(tsp_instance* tsp_in)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		tsp_in->sol[i] = i;
	}
	tsp_in->sol[tsp_in->num_nodes] = tsp_in->sol[0];

	evaluate_sol(tsp_in);
	print_cost(tsp_in);

	if (tsp_in->plot)
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