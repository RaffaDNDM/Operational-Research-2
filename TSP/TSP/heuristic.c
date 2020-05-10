#include "heuristic.h"
#include "utility.h"

void heuristic_solver(tsp_instance* tsp_in)
{
	switch (tsp_in->alg)
	{
		case 6:
		{
			printf(LINE);
			printf("Heuristic nearest neighborhood solver\n");
			time_t start = clock();
			nearest_neighborhood(tsp_in);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}
		case 7:
		{

			printf(LINE);
			printf("Heuristic insertion solver\n");
			time_t start = clock();
			insertion(tsp_in);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}
	}

	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = calloc(tsp_in->num_nodes, sizeof(int));

	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
			case 6: case 7:
			{
				define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
				break;
			}
		}

		plot(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
	}

	free(tsp_in->sol);
}

void nearest_neighborhood(tsp_instance* tsp_in)
{
	int i = 0;

	int* nodes =(int *) calloc(tsp_in->num_nodes, sizeof(int));

	int num_edges = tsp_in->num_nodes * (tsp_in->num_nodes - 1) / 2;
	tsp_in->sol = (double*)calloc(((size_t)num_edges), sizeof(double));

	if (tsp_in->integerDist)
		tsp_in->bestCostI = 0;
	else
		tsp_in->bestCostD = 0.0;

	nodes[0] = 1;

	for (; i < (tsp_in->num_nodes); i++)
	{
		double min_dist = DBL_MAX;
		int best = tsp_in->num_nodes;

		min_cost(tsp_in, nodes, i, &min_dist, &best);
		
		if (best == tsp_in->num_nodes)
		{
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, 0, tsp_in, &x);
				min_dist = (double)x;
			}
			else
				dist(i, 0, tsp_in, &min_dist);

			(tsp_in->sol)[xpos(tsp_in, i, 0)] = 1.0;
		}
		else
			(tsp_in->sol)[xpos(tsp_in, i, best)] = 1.0;

		if (tsp_in->integerDist)
			tsp_in->bestCostI += (int) min_dist;
		else
			tsp_in->bestCostD += min_dist;

		i = best-1;
	}

}

void insertion(tsp_instance* tsp_in)
{
	int num_edges = tsp_in->num_nodes * (tsp_in->num_nodes-1) /2 ;
	tsp_in->sol = (double*)calloc(((size_t)num_edges), sizeof(double));

	double max_dist = -1;
	int indices[2];
	int i = 0;

	for (; i < tsp_in->num_nodes; i++)
	{
		int j = i + 1;
		for (; j < tsp_in->num_nodes; j++)
		{
			double c;

			if (tsp_in->integerDist)
			{
				int x;
				dist(i, j, tsp_in, &x);
				c = (double)x;
			}
			else
				dist(i, j, tsp_in, &c);

			if (c > max_dist)
			{
				max_dist = c;
				indices[0] = i;
				indices[1] = j;
			}
		}
	}

	int* node1 = (int *) calloc((size_t) tsp_in->num_nodes, sizeof(int));
	int* node2 = (int *) calloc((size_t) tsp_in->num_nodes, sizeof(int));
	int* visited_nodes = (int*) calloc((size_t) tsp_in->num_nodes, sizeof(int));
	double* costs = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));

	int count = 2;

	visited_nodes[0] = indices[0];
	visited_nodes[1] = indices[1];

	node1[0] = indices[0];
	node2[0] = indices[1];
	node1[1] = indices[0];
	node2[1] = indices[1];

	costs[0] = costs[1] = max_dist;

	double best_cost = max_dist * 2;
	int i_best;
	int k_best;
	

	for (; count < tsp_in->num_nodes - 1; count ++)
	{
		double best_cost_h = DBL_MAX;

		min_extra_mileage(tsp_in, count, visited_nodes, node1, node2, costs, &i_best, &k_best ,&best_cost_h, &best_cost);
		

		if (tsp_in->integerDist)
		{
			int x1, x2;
			dist(node1[k_best], node2[k_best], tsp_in, &x1);
			dist(node1[count], node2[count], tsp_in, &x2);
			costs[k_best] = (double) x1;
			costs[count] = (double) x2;
		}
		else
		{
			dist(node1[k_best], node2[k_best], tsp_in, &(costs[k_best]));
			dist(node1[count], node2[count], tsp_in, &(costs[count]));
		}
	}
	
	for (i = 0; i < tsp_in->num_nodes-1; i++)
		(tsp_in->sol)[xpos(tsp_in, node1[i], node2[i])] = 1.0;

	if (tsp_in->integerDist)
		tsp_in->bestCostI = (int)(best_cost + CAST_PRECISION);
	else
		tsp_in->bestCostD = best_cost;
}

void min_cost(tsp_instance* tsp_in, int* nodes, int i, double* min_dist, int* best)
{
#ifdef GRASP

	srand(time(NULL));
	double min[] = { DBL_MAX , DBL_MAX, DBL_MAX };
	int min_pos[3];
	min_pos[0] = tsp_in->num_nodes;
	min_pos[1] = tsp_in->num_nodes;
	min_pos[2] = tsp_in->num_nodes;
	int j = 0;
	for (; j < tsp_in->num_nodes; j++)
	{
		if (i == j)
			continue;

		if (nodes[j])
			continue;

		double c; //cost of the edge

		if (tsp_in->integerDist)
		{
			int x;
			dist(i, j, tsp_in, &x);
			c = (double)x;
		}
		else
			dist(i, j, tsp_in, &c);

		if (c < min[0])
		{
			min[2] = min[1];
			min[1] = min[0];
			min[0] = c;
			min_pos[0] = j;
		}
		else if (c < min[1])
		{
			min[2] = min[1];
			min[1] = c;
			min_pos[1] = j;
		}
		else if (c < min[2])
		{
			min[2] = c;
			min_pos[2] = j;
		}

	}

	int n = rand() % 8;
	if (n < 3)
	{
		*min_dist = min[0];
		*best = min_pos[0];
	}
	else if (n < 6)
	{
		*min_dist = min[1];
		*best = min_pos[1];
	}
	else
	{
		*min_dist = min[2];
		*best = min_pos[2];
	}

	#else

	int j = 0;
	for (; j < tsp_in->num_nodes; j++)
	{
		if (i == j)
			continue;

		if (nodes[j])
			continue;

		double c; //cost of the edge

		if (tsp_in->integerDist)
		{
			int x;
			dist(i, j, tsp_in, &x);
			c = (double)x;
		}
		else
			dist(i, j, tsp_in, &c);

		if (c < (*min_dist) )
		{
			(*min_dist) = c;
			(*best) = j;
		}
	}
	#endif

	nodes[*best] = 1;

}

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost )
{
#ifdef GRASP

	printf("count %d\n", count);

	srand(time(NULL));
	double min[] = { DBL_MAX , DBL_MAX, DBL_MAX };
	int min_nodes[3];
	int min_edges[3];
	int h = 0;
	for (; h < tsp_in->num_nodes; h++)
	{
		int k = 0;
		int jump = 0;
		for (; k < count; k++)
		{
			if (h == visited_nodes[k])
			{
				jump = 1;
				break;
			}
		}

		double min_h = DBL_MAX;
		int k_h;

		for (k = 0; k < count && jump == 0; k++)
		{
			double c;

			if (tsp_in->integerDist)
			{
				int x1, x2;
				dist(h, node1[k], tsp_in, &x1);
				dist(h, node2[k], tsp_in, &x2);
				c = (double)(x1 + x2) - costs[k];
			}
			else
			{
				double x1, x2;
				dist(h, node1[k], tsp_in, &x1);
				dist(h, node2[k], tsp_in, &x2);
				c = x1 + x2 - costs[k];
			}

			if (c < min_h)
			{
				min_h = c;
				k_h = k;
			}
		}

		if (min_h < min[0])
		{
			min[2] = min[1];
			min[1] = min[0];
			min[0] = min_h;

			min_edges[2] = min_edges[1];
			min_edges[1] = min_edges[0];
			min_edges[0] = k_h;

			min_nodes[2] = min_nodes[1];
			min_nodes[1] = min_nodes[0];
			min_nodes[0] = h;
		}
		else if (min_h < min[1])
		{
			min[2] = min[1];
			min[1] = min_h;

			min_edges[2] = min_edges[1];
			min_edges[1] = k_h;

			min_nodes[2] = min_nodes[1];
			min_nodes[1] = h;
		}
		else if (min_h < min[2])
		{
			min[2] = min_h;

			min_edges[2] = k_h;

			min_nodes[2] = h;
		}
	}
	
	
	if (count < tsp_in->num_nodes - 3)
	{
		int n = rand() % 9;
		if (n < 3)
		{
			(*best_cost_h) = min[0];
			(*k_best) = min_edges[0];
			(*i_best) = min_nodes[0];
		}
		else if (n < 6)
		{
			(*best_cost_h) = min[1];
			(*k_best) = min_edges[1];
			(*i_best) = min_nodes[1];
		}
		else
		{
			(*best_cost_h) = min[2];
			(*k_best) = min_edges[2];
			(*i_best) = min_nodes[2];
		}
	}
	else
	{
		(*best_cost_h) = min[0];
		(*i_best) = min_nodes[0];
		(*k_best) = min_edges[0];
	}

#else
	
	int h = 0;
	for (; h < tsp_in->num_nodes; h++)
	{
		int k = 0;
		int jump = 0;
		for (; k < count; k++)
		{
			if (h == visited_nodes[k])
			{
				jump = 1;
				break;
			}
		}

		double min_h = DBL_MAX;
		int k_h;

		for (k = 0; k < count && jump == 0; k++)
		{
			double c;

			if (tsp_in->integerDist)
			{
				int x1, x2;
				dist(h, node1[k], tsp_in, &x1);
				dist(h, node2[k], tsp_in, &x2);
				c = (double)(x1 + x2) - costs[k];
			}
			else
			{
				double x1, x2;
				dist(h, node1[k], tsp_in, &x1);
				dist(h, node2[k], tsp_in, &x2);
				c = x1 + x2 - costs[k];
			}

			if (c < min_h)
			{
				min_h = c;
				k_h = k;
			}
		}

		if (min_h < (*best_cost_h) )
		{
			(*best_cost_h) = min_h;
			(*i_best) = h;
			(*k_best) = k_h;
		}
	}
#endif

	(*best_cost) += (*best_cost_h);

	visited_nodes[count] = (*i_best);

	int temp = node2[(*k_best)];
	node2[(*k_best)] = (*i_best);
	node1[count] = (*i_best);
	node2[count] = temp;
}