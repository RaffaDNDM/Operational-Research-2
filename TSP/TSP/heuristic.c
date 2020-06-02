/**
	@file heuristic.c
	@author Cristina Fabris
	@author Raffaele Di Nardo Di Maio
	@brief Heuristic solvers.
*/

#include "heuristic.h"
#include "utility.h"

pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;


void heuristic_solver(tsp_instance* tsp_in)
{
	time_t start = clock();

	tsp_in->bestCostD = DBL_MAX;
	tsp_in->bestCostI = INT_MAX;

	int  num_edges = (tsp_in->num_nodes * (tsp_in->num_nodes)) / 2;
	tsp_in->sol = (double*) calloc((size_t) num_edges, sizeof(double));
	int* succ = (int*) calloc((size_t)tsp_in->num_nodes, sizeof(int));

	printf(LINE);
	printf("[Construction] ");

	if (!CONSTRUCTION_TYPE)
	{
		printf("Nearest Neighborhood\n");
	}
	else
	{
		printf("Insertion\n");
	}


	printf("[Meta-heuristic] ");

	switch (tsp_in->alg)
	{
		case 7:
		{
			printf("VNS\n");
			break;
		}
		case 8:
		{
			printf("Tabu search\n");
			break;
		}
	}
	printf(LINE);

	#ifdef MULTI_START
		pthread_t threads[NUM_MULTI_START];
		thread_args param[NUM_MULTI_START];

		int i = 0;
		for (; i < NUM_MULTI_START; i++)
		{
			param[i].tsp_in = tsp_in;
			param[i].succ = succ;
			param[i].seed = STEP_SEED*(i+1);

			pthread_create(&threads[i], NULL, computeSolution, (void*)&param[i]);
		}

		for (i = 0; i < NUM_MULTI_START; i++)
		{
			int rc = pthread_join(threads[i], NULL);

			if (rc)
				exit(-1);
		}

	#else
		pthread_t thread;
		thread_args param;

		param.tsp_in = tsp_in;
		param.succ = succ;
		param.seed = 0;

		pthread_create(&thread, NULL, computeSolution, (void*)&param);

		int rc = pthread_join(thread, NULL);

		if (rc)
			exit(-1);

	#endif

	time_t end = clock();
	tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	int* comp = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
		case 7: case 8:
		{
			int k;
			for (k = 0; k < tsp_in->num_nodes; k++)
				comp[k] = 1;
			break;
		}
		}

		plot(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
	}

	free(tsp_in->sol);
}

void* computeSolution(void* param)
{
	thread_args* args = (thread_args*) param;
	double best_cost = 0.0;
	int* visited_nodes = (int*)calloc((size_t) args->tsp_in->num_nodes, sizeof(int));

	if (!CONSTRUCTION_TYPE)
		nearest_neighborhood(args->tsp_in, visited_nodes, &best_cost, args->seed);
	else
		insertion(args->tsp_in, visited_nodes, &best_cost, args->seed);

	switch (args->tsp_in->alg)
	{
		case 7:
		{
			vns(args->tsp_in, visited_nodes, &best_cost);
			break;
		}
		case 8:
		{
			tabu_search(args->tsp_in, visited_nodes, &best_cost);
			break;
		}
	}

	pthread_mutex_lock(&mutex);
	printf("Cost: %lf\n", best_cost);

	if (args->tsp_in->integerDist)
	{
		if (((int)best_cost) < args->tsp_in->bestCostI)
		{
			args->tsp_in->bestCostI = (int)best_cost;
			update_solution(visited_nodes, args->tsp_in->sol, args->tsp_in->num_nodes);
			succ_construction(visited_nodes, args->succ, args->tsp_in->num_nodes);
		}
	}
	else
	{
		if (best_cost < args->tsp_in->bestCostD)
		{
			args->tsp_in->bestCostD = best_cost;
			update_solution(visited_nodes, args->tsp_in->sol, args->tsp_in->num_nodes);
			succ_construction(visited_nodes, args->succ, args->tsp_in->num_nodes);
		}
	}

	pthread_mutex_unlock(&mutex);

	free(visited_nodes);
	pthread_exit(NULL);
}

void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, int seed)
{
	int i = 0;

	int* nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	//vettore di dim num nodes -> nodes utilizzato per fare le verifiche in min_cost
	//metto a 1 i nodi visitati

	(*best_cost) = 0.0;

	visited_nodes[0] = 0;
	nodes[0] = 1;
	int count = 1;

	for (; i < (tsp_in->num_nodes); i++)
	{
		double min_dist = DBL_MAX;
		int best = tsp_in->num_nodes;

		min_cost(tsp_in, nodes, i, &min_dist, &best, seed);

		if (best == tsp_in->num_nodes)
		{
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, 0, tsp_in, &x);
				min_dist = (double) x;
			}
			else
				dist(i, 0, tsp_in, &min_dist);

			//printf("count = %d\n", count);

		}
		else
		{
			visited_nodes[count] = best;//aggiungere nodo a visited nodes
			count++;
		}

		(*best_cost) += min_dist;

		i = best-1;
	}

}

void insertion(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, int seed)
{
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
	double* costs = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));

	int count = 2;

	visited_nodes[0] = indices[0];
	visited_nodes[1] = indices[1];

	node1[0] = indices[0];
	node2[0] = indices[1];
	node1[1] = indices[1];
	node2[1] = indices[0];

	costs[0] = costs[1] = max_dist;

	(*best_cost) = max_dist * 2;

	int i_best;

	for (; count < tsp_in->num_nodes ; count++)
	{
		double best_cost_h = DBL_MAX;
		int k_best;

		min_extra_mileage(tsp_in, count, visited_nodes, node1, node2, costs, &i_best, &k_best ,&best_cost_h, best_cost, seed);

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
		/*
		int* succ_nodes = calloc((size_t)tsp_in->num_nodes, sizeof(int));
		succ_construction(visited_nodes, succ_nodes, tsp_in->num_nodes);
		update_solution(visited_nodes, tsp_in, tsp_in->num_nodes);
		int* comp = calloc((size_t)tsp_in->num_nodes, sizeof(int));
		int s;
		for (s = 0; s < tsp_in->num_nodes; s++)
			comp[s] = 1;
		int n_comps = 1;
		plot(tsp_in, succ_nodes, comp, &n_comps);
		*/
	}
}

void min_cost(tsp_instance* tsp_in, int* nodes, int i, double* min_dist, int* best, int seed)
{
#ifdef GRASP

	if(!seed)
		srand(time(NULL));
	else
		srand(seed);

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

			min_pos[2] = min_pos[1];
			min_pos[1] = min_pos[0];
			min_pos[0] = j;
		}
		else if (c < min[1])
		{
			min[2] = min[1];
			min[1] = c;

			min_pos[2] = min_pos[1];
			min_pos[1] = j;
		}
		else if (c < min[2])
		{
			min[2] = c;
			min_pos[2] = j;
		}

	}


	int max = 0;
	if (min[1] == DBL_MAX)
		max = 3;
	else if (min[2] == DBL_MAX)
		max = 6;
	else
		max = 9;

	int n = rand() % max;

	if (n < 3 )
	{
		*min_dist = min[0];
		*best = min_pos[0];
	}
	else if (n < 6 )
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

	if(*best != tsp_in->num_nodes)
		nodes[*best] = 1;

}

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost, int seed)
{
#ifdef GRASP

	if (!seed)
		srand(time(NULL));
	else
		srand(seed);

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

		for (k = 0; k < count && !jump; k++)
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


	if (count < tsp_in->num_nodes - 2)
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

		for (k = 0; k < count && !jump; k++)
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

	//visited_nodes[count] = (*i_best);
	int j = 0;
	while (visited_nodes[j] != node1[*k_best])
		j++;

	int s = count;
	for (; s > j + 1; s--)
		visited_nodes[s] = visited_nodes[s - 1];

	visited_nodes[s] = *i_best;

	int temp = node2[(*k_best)];
	node2[(*k_best)] = (*i_best);
	node1[count] = (*i_best);
	node2[count] = temp;
	/*
	printf("visited notes: ");
	s = 0;
	for (; s <= count; s++)
		printf(" %d ", visited_nodes[s]);

	printf("\nnode1: ");
	for (s = 0; s <= count; s++)
		printf(" %d ", node1[s]);


	printf("\nnode2: ");
	for (s = 0; s <= count; s++)
		printf(" %d ", node2[s]);

	printf("\n%s", LINE);

	s = 2;*/
}

void greedy_refinement(tsp_instance* tsp_in, int* visited_nodes, double* best_cost)
{
	static int greedy_count = 0;
	greedy_count++;
	int* succ = calloc((size_t)tsp_in->num_nodes, sizeof(int));
	succ_construction(visited_nodes, succ, tsp_in->num_nodes);

	double check_cost;

	do
	{
		check_cost = (*best_cost);

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)
		{
			double cost_i_k; //cost [i, succ[i]]
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, succ[i], tsp_in, &x);
				cost_i_k = (double)(x);
			}
			else
			{
				dist(i, succ[i], tsp_in, &cost_i_k);
			}

			int j = 0;
			for (; j < tsp_in->num_nodes; j++)
			{
				if (j != i && j != succ[i] && succ[j] != i && succ[j] != succ[i])
				{
					double cost_j_h; //cost[i, succ[i]]
					double cost_i_j;
					double cost_k_h;

					if (tsp_in->integerDist)
					{
						int x1, x2, x3;
						dist(j, succ[j], tsp_in, &x1);
						dist(i, j, tsp_in, &x2);
						dist(succ[i], succ[j], tsp_in, &x3);
						cost_j_h = (double)x1;
						cost_i_j = (double)x2;
						cost_k_h = (double)x3;
					}
					else
					{
						dist(j, succ[j], tsp_in, &cost_j_h);
						dist(i, j, tsp_in, &cost_i_j);
						dist(succ[i], succ[j], tsp_in, &cost_k_h);
					}

					double delta = cost_i_j + cost_k_h - cost_i_k - cost_j_h;

					if (delta < 0.0)
					{
						(*best_cost) += delta;

						int k = succ[i];
						int next = succ[i];

						int count = 0;
						int* orientation = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

						while (next != j)
						{
							orientation[count++] = next;
							next = succ[next];
						}

						next = succ[k];
						int next_count = 0;

						while (next != j)
						{
							int temp = succ[next];
							succ[next] = orientation[next_count++];
							next = temp;
						}

						succ[i] = j;
						succ[k] = succ[j];
						succ[j] = orientation[count - 1];

						free(orientation);
						break;
					}
				}
			}
		}
	}
	while (check_cost != (*best_cost)); //abs(check_cost - (*cost))<1e-10

	visited_nodes[0] = 0;
	int i = 0;
	int count = 1;
	int node = -1;

	for(; count<tsp_in->num_nodes; count++)
	{
		node = succ[i];
		i = node;
		visited_nodes[count] = node;  //inserirlo nell'algoritmo
		//printf("[ %d ] count: %d\n", greedy_count, visited_nodes[count]);
	}

	free(succ);
}

void vns(tsp_instance* tsp_in, int* visited_nodes, double *best_cost)
{
	srand(time(NULL));

	int* local_min_visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
		local_min_visited_nodes[i] = visited_nodes[i];

	int num_local_mins = 0;
	int num_iterations = 0;
	double local_min_cost = (*best_cost);


	while(num_local_mins < MAX_LOCAL_MINS && num_iterations < MAX_NUM_ITERATIONS)
	{
		int* kopt_visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

		for (i=0; i < tsp_in->num_nodes; i++)
		{
			kopt_visited_nodes[i] = local_min_visited_nodes[i];
		}

		/*
		for (i = 0; i < tsp_in->num_nodes; i++)
			printf("[%d]  %d\n", i, kopt_visited_nodes[i]);
		*/

		int t = rand() % tsp_in->num_nodes;
		int k = 1;
		int max_k = tsp_in->num_nodes / 2;
		int found_min = 0;
		double kopt_cost = local_min_cost;

		for (; k < max_k && !found_min; k++)
		{
			num_iterations++;
			int next_index = (t + k) % tsp_in->num_nodes;
			int temp = kopt_visited_nodes[next_index];
			kopt_visited_nodes[next_index] = kopt_visited_nodes[t];
			kopt_visited_nodes[t] = temp;

			if (tsp_in->integerDist)
			{
				int c_old[4], c_new[4];

				dist(kopt_visited_nodes[(t - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], kopt_visited_nodes[next_index], tsp_in, &c_old[0]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(t + 1) % tsp_in->num_nodes], tsp_in, &c_old[1]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(next_index + 1) % tsp_in->num_nodes], tsp_in, &c_old[2]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(next_index - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], tsp_in, &c_old[3]);
				dist(kopt_visited_nodes[(t - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], kopt_visited_nodes[t], tsp_in, &c_new[0]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(t + 1) % tsp_in->num_nodes], tsp_in, &c_new[1]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(next_index + 1) % tsp_in->num_nodes], tsp_in, &c_new[2]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(next_index - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], tsp_in, &c_new[3]);

				for (i = 0; i < 4; i++)
					kopt_cost += (double) (c_new[i]-c_old[i]);
			}
			else
			{
				double c_old[4], c_new[4];

				dist(kopt_visited_nodes[(t - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], kopt_visited_nodes[next_index], tsp_in, &c_old[0]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(t + 1) % tsp_in->num_nodes], tsp_in, &c_old[1]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(next_index + 1) % tsp_in->num_nodes], tsp_in, &c_old[2]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(next_index - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], tsp_in, &c_old[3]);
				dist(kopt_visited_nodes[(t - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], kopt_visited_nodes[t], tsp_in, &c_new[0]);
				dist(kopt_visited_nodes[t], kopt_visited_nodes[(t + 1) % tsp_in->num_nodes], tsp_in, &c_new[1]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(next_index + 1) % tsp_in->num_nodes], tsp_in, &c_new[2]);
				dist(kopt_visited_nodes[next_index], kopt_visited_nodes[(next_index - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], tsp_in, &c_new[3]);

				for (i = 0; i < 4; i++)
					kopt_cost += (c_new[i] - c_old[i]);
			}


			int* result_kopt_visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));


			for(i = 0; i < tsp_in->num_nodes; i++)
			{
				result_kopt_visited_nodes[i] = kopt_visited_nodes[i];
			}

			double result_kopt_cost = kopt_cost;
			greedy_refinement(tsp_in, result_kopt_visited_nodes, &result_kopt_cost);

			if(local_min_cost != result_kopt_cost)//abs(local_min_cost - result_kopt_cost)>1e-10
			{
				num_local_mins++;
				local_min_cost = result_kopt_cost;

				for (i = 0; i < tsp_in->num_nodes; i++)
				{
					local_min_visited_nodes[i] = result_kopt_visited_nodes[i];
				}

				if (local_min_cost < (*best_cost))
				{
					found_min = 1;
					(*best_cost) = local_min_cost;
					for (i = 0; i < tsp_in->num_nodes; i++)
					{
						visited_nodes[i] = result_kopt_visited_nodes[i];
					}
				}
			}

			free(result_kopt_visited_nodes);

			if (num_iterations == MAX_NUM_ITERATIONS || num_local_mins == MAX_LOCAL_MINS)
				break;
		}

		free(kopt_visited_nodes);
	}
}

void update_solution(int* visited_nodes, double* sol, int num_nodes)//necessario aver gi� allocato sol , un vettore di double generico
{
	int num_edges = num_nodes * (num_nodes - 1) / 2;

	int i;
	for (i = 0; i < num_edges; i++)
		sol[i] = 0.0;

	for (i = 0; i < num_nodes; i++)
		sol[generic_xpos(visited_nodes[i], visited_nodes[(i + 1) % num_nodes], num_nodes)] = 1.0;

}

void succ_construction(int* visited_nodes, int* succ, int num_nodes)//succ deve essere gi� allocato
{
	int i;
	for (i = 0; i < num_nodes; i++)
		succ[visited_nodes[i]] = visited_nodes[(i+1) % num_nodes];

}

void tabu_search(tsp_instance* tsp_in, int* visited_nodes, double* best_cost)
{
	greedy_refinement(tsp_in, visited_nodes, best_cost); //prima di iniziare l'algoritmo devo essere sicura di trovarmi in un minimo locale

	int** tabu_list = (int**)calloc((size_t)2, sizeof(int*));//[tabu_list[0][i], tabu_list[1][i]]
	if (tabu_list == NULL)
		printf("ERROR\n");

	int tenure = ceil(tsp_in->num_nodes / 10);


	tabu_list[0] = (int*)calloc((size_t)tenure, sizeof(int));
	if (tabu_list[0] == NULL)
		printf("ERROR\n");
	tabu_list[1] = (int*)calloc((size_t)tenure, sizeof(int));
	if (tabu_list[1] == NULL)
		printf("ERROR\n");

	int num_tabu_edges = 0;

	//in position i has the cost of the edge [visited_nodes[i], visited_nodes[i+1] ]
	double* costs = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));
	if (costs == NULL)
		printf("ERROR\n");
	int i;
	for (i = 0; i < tsp_in->num_nodes - 1; i++)
	{
		if (tsp_in->integerDist)
		{
			int x;
			dist(visited_nodes[i], visited_nodes[i + 1], tsp_in, &x);
			costs[i] = (double)x;
		}
		else
			dist(visited_nodes[i], visited_nodes[i + 1], tsp_in, &(costs[i]));
	}

	int* nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	if (nodes == NULL)
		printf("ERROR\n");
	int j;
	for (j = 0; j < tsp_in->num_nodes; j++)
		nodes[j] = visited_nodes[j];

	double actual_cost = *best_cost;

	int num_iteration = 0;
	for (; num_iteration < MAX_NUM_ITERATIONS; num_iteration++)
	{
		double min_increase = DBL_MAX;
		int pos1 = -1;
		int pos2 = -1;
		int restart = 1;
		int refinement_done = 0;
		int negative = 0;

		for (i = 0; i < tsp_in->num_nodes && !negative; i++)
		{
			int start_edge1, end_edge1, start_edge2, end_edge2;

			for (j = 1; j < 3; j++)
			{
				int jump = 0;

				start_edge1 = nodes[(i - 1 + tsp_in->num_nodes) % tsp_in->num_nodes];
				end_edge1 = nodes[(i + j) % tsp_in->num_nodes];
				start_edge2 = nodes[i];
				end_edge2 = nodes[(i + j + 1) % tsp_in->num_nodes];

				int k;
				for (k = 0; k < num_tabu_edges; k++)
				{
					if ((start_edge1 == tabu_list[0][k] && end_edge1 == tabu_list[1][k])
						|| (start_edge2 == tabu_list[0][k] && end_edge2 == tabu_list[1][k])
						|| (start_edge1 == tabu_list[1][k] && end_edge1 == tabu_list[0][k])
						|| (start_edge2 == tabu_list[1][k] && end_edge2 == tabu_list[0][k])
						)
					{
						jump = 1;
						break;
					}
				}

				if (jump)
					continue;

				double delta;
				if (tsp_in->integerDist)
				{
					int x1, x2;
					dist(start_edge1, end_edge1, tsp_in, &x1);
					dist(start_edge2, end_edge2, tsp_in, &x2);
					delta = (double)x1 + x2 - costs[(i + j) % tsp_in->num_nodes] - costs[i - 1];
				}
				else
				{
					double x1, x2;
					dist(start_edge1, end_edge1, tsp_in, &x1);
					dist(start_edge2, end_edge2, tsp_in, &x2);
					delta = x1 + x2 - costs[(i + j) % tsp_in->num_nodes] - costs[i - 1];
				}

				if (abs(delta) < 1e-10)//in valore assoluto
					delta = 0.0;

				if (delta < 0)
				{
					negative = 1;
					break;
				}

				if (delta < min_increase && delta != 0)
					//if (min_increase - delta > EPS && delta != 0.0)
				{
					min_increase = delta;
					pos1 = i;
					pos2 = (i + j) % tsp_in->num_nodes;
				}
			}

		}

		//if (num_negative_delta == 2 * tsp_in->num_nodes)
		if (negative)
		{
			greedy_refinement_for_tabu_search(tsp_in, nodes, tabu_list, tenure, &num_tabu_edges, &actual_cost);

		}
		else
		{
			int temp = nodes[pos1];
			nodes[pos1] = nodes[pos2];
			nodes[pos2] = temp;

			add_element(tabu_list[0], tabu_list[1], tenure, nodes[(pos1 - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], nodes[pos2]);
			add_element(tabu_list[0], tabu_list[1], tenure, nodes[pos1], nodes[(pos2 + 1) % tsp_in->num_nodes]);

			if (num_tabu_edges < tenure)
				num_tabu_edges = num_tabu_edges + 2;

			actual_cost = actual_cost + min_increase;
		}

		//aggionare vettore di costi ****************FARLO IN MANIERA PI� DECENTE***************
		int h;
		for (h = 0; h < tsp_in->num_nodes - 1; h++)
		{
			if (tsp_in->integerDist)
			{
				int x;
				dist(nodes[h], nodes[h + 1], tsp_in, &x);
				costs[h] = (double)x;
			}
			else
				dist(nodes[h], nodes[h + 1], tsp_in, &(costs[h]));
		}

		if (actual_cost < *best_cost)//controllare precisione
		{
			*best_cost = actual_cost;

			for (j = 0; j < tsp_in->num_nodes; j++)
				visited_nodes[j] = nodes[j];
		}
	}

	free(costs);
	free(tabu_list[0]);
	free(tabu_list[1]);
	free(tabu_list);
	free(nodes);
}

static int end_list = 0;

void add_element(int* list1, int* list2, int dimension, int element1, int element2)
{
	list1[end_list] = element1;
	list2[end_list] = element2;
	end_list = (end_list + 1) % dimension;
}

void greedy_refinement_for_tabu_search(tsp_instance* tsp_in, int* visited_nodes, int** tabu_list, int tenure, int* num_tabu_edges, double* best_cost)
{
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));
	succ_construction(visited_nodes, succ, tsp_in->num_nodes);

	double check_cost;

	do
	{
		check_cost = (*best_cost);

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)
		{
			double cost_i_k; //cost [i, succ[i]]
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, succ[i], tsp_in, &x);
				cost_i_k = (double)(x);
			}
			else
			{
				dist(i, succ[i], tsp_in, &cost_i_k);
			}

			int j = 0;
			for (; j < tsp_in->num_nodes; j++)
			{
				if (j != i && j != succ[i] && succ[j] != i && succ[j] != succ[i])
				{
					int k;
					int is_tabu = 0;
					for (k = 0; k < *num_tabu_edges; k++)
					{
						if ((i == tabu_list[0][k] && j == tabu_list[1][k])
							|| (succ[i] == tabu_list[0][k] && succ[j] == tabu_list[1][k])
							|| (i == tabu_list[1][k] && j == tabu_list[0][k])
							|| (succ[i] == tabu_list[1][k] && succ[j] == tabu_list[0][k]))
						{
							is_tabu = 1;
							break;
						}
					}

					if (!is_tabu)
					{
						double cost_j_h; //cost[i, succ[i]]
						double cost_i_j; //cost[i,j]
						double cost_k_h; //cost[succ[i], succ[j]]

						if (tsp_in->integerDist)
						{
							int x1, x2, x3;
							dist(j, succ[j], tsp_in, &x1);
							dist(i, j, tsp_in, &x2);
							dist(succ[i], succ[j], tsp_in, &x3);
							cost_j_h = (double)x1;
							cost_i_j = (double)x2;
							cost_k_h = (double)x3;
						}
						else
						{
							dist(j, succ[j], tsp_in, &cost_j_h);
							dist(i, j, tsp_in, &cost_i_j);
							dist(succ[i], succ[j], tsp_in, &cost_k_h);
						}

						double delta = cost_i_j + cost_k_h - cost_i_k - cost_j_h;

						if (delta < 0.0) //if (0.0 - delta > EPS)
						{
							add_element(tabu_list[0], tabu_list[1], tenure, i, succ[i]);
							add_element(tabu_list[0], tabu_list[1], tenure, j, succ[j]);

							if (*num_tabu_edges < tenure)
								*num_tabu_edges = *num_tabu_edges + 2;

							(*best_cost) += delta;

							int k = succ[i];
							int next = succ[i];

							int count = 0;
							int* orientation = (int*)calloc((size_t)tsp_in->num_nodes,sizeof(int));

							while (next != j)
							{
								orientation[count++] = next;
								next = succ[next];
							}

							next = succ[k];
							int next_count = 0;

							while (next != j)
							{
								int temp = succ[next];
								succ[next] = orientation[next_count++];
								next = temp;
							}

							succ[i] = j;
							succ[k] = succ[j];
							succ[j] = orientation[count - 1];

							free(orientation);
							break;
						}
					}
				}
			}
		}

	} while (abs(check_cost - (*best_cost)) > 1e-10);


	visited_nodes[0] = 0;
	int i = 0;
	int count = 1;
	int node = -1;

	for (; count < tsp_in->num_nodes; count++)
	{
		node = succ[i];
		i = node;
		visited_nodes[count] = node;  //inserirlo nell'algoritmo
	}

	free(succ);
}
