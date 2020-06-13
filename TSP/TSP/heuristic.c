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

	int  num_edges = (tsp_in->num_nodes) * (tsp_in->num_nodes-1) / 2;
	tsp_in->sol = (double*) calloc((size_t) num_edges, sizeof(double));
	int* succ = (int*) calloc((size_t)tsp_in->num_nodes, sizeof(int));

	printf("%sHeuristic solver%s\n",RED, WHITE);
	printf("%s[Construction]%s ", BLUE, WHITE);

	if (!CONSTRUCTION_TYPE)
	{
		printf("Nearest Neighborhood\n");
	}
	else
	{
		printf("Insertion\n");
	}


	printf("%s[Meta-heuristic]%s ", BLUE, WHITE);

	switch (tsp_in->alg)
	{
		case 8:
		{
			printf("VNS\n");
			break;
		}
		case 9:
		{
			printf("Tabu search\n");
			break;
		}
		case 10:
		{
			printf("Simulated anealling\n");
			break;
		}
		case 11:
		{
			printf("Genetic\n");
			break;
		}
	}
	printf("%s%s%s", RED, LINE, WHITE);

	#ifdef MULTI_START
		pthread_t threads[NUM_MULTI_START];
		thread_args param[NUM_MULTI_START];

		int i = 0;
		for (; i < NUM_MULTI_START; i++)
		{
			param[i].tsp_in = tsp_in;
			param[i].succ = succ;
			param[i].seed = STEP_SEED*(i+1);
			param[i].start = start;

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
		param.start = start;

		pthread_create(&thread, NULL, computeSolution, (void*)&param);

		int rc = pthread_join(thread, NULL);

		if (rc)
			exit(-1);

	#endif

	time_t end = clock();
	tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
	print_cost(tsp_in);
	printf("%sExecution time:%s %.3lf seconds\n",GREEN, WHITE, tsp_in->execution_time);
	printf("%s%s%s", RED, LINE, WHITE);

	int* comp = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
		case 8: case 9: case 10:
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

	time_t start = clock();
	double remaining_time = args->tsp_in->deadline -( (double)(start - args->start) / (double) CLOCKS_PER_SEC );

	if (!CONSTRUCTION_TYPE)
		nearest_neighborhood(args->tsp_in, visited_nodes, &best_cost, args->seed, -1);
	else
		insertion(args->tsp_in, visited_nodes, &best_cost, args->seed);

	greedy_refinement(args->tsp_in, visited_nodes, &best_cost);

	time_t end_first_generation = clock();
	remaining_time = remaining_time - ((double)(end_first_generation-start)/(double)CLOCKS_PER_SEC);

	switch (args->tsp_in->alg)
	{
		case 8:
		{
			vns(args->tsp_in, visited_nodes, &best_cost, remaining_time);
			break;
		}
		case 9:
		{
			tabu_search(args->tsp_in, visited_nodes, &best_cost, remaining_time);
			break;
		}
		case 10:
		{
			simulated_annealing(args->tsp_in, visited_nodes, &best_cost, remaining_time);
			break;
		}
	}

	pthread_mutex_lock(&mutex);
	printf("%sCost:%s %.2lf\n",GREEN, WHITE, best_cost);

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

void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, int seed, int first_node) 
{
	int* nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	//vettore di dim num nodes -> nodes utilizzato per fare le verifiche in min_cost
	//metto a 1 i nodi visitati

	(*best_cost) = 0.0;

	if (tsp_in->alg == 11)
	{
		visited_nodes[0] = first_node;
		nodes[visited_nodes[0]] = 1;
		//i = first_node;//aggiunto
	}
	else
	{
		visited_nodes[0] = 0;
		nodes[0] = 1;
	}

	int count = 1;

	int i = visited_nodes[0];
	int iter = 0;
	for (; i < (tsp_in->num_nodes); i++)
	//for(; iter< tsp_in->num_nodes; iter++)
	{
		double min_dist = DBL_MAX;
		int best = tsp_in->num_nodes;

		
		min_cost(tsp_in, nodes, i, &min_dist, &best, seed);

		if (best == tsp_in->num_nodes)
		{
			if (tsp_in->integerDist)
			{
				int x;
				dist(i, visited_nodes[0], tsp_in, &x);
				min_dist = (double) x;
			}
			else
				dist(i, visited_nodes[0], tsp_in, &min_dist);

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
	int* succ = calloc((size_t)tsp_in->num_nodes, sizeof(int));
	succ_construction(visited_nodes, succ, tsp_in->num_nodes);
	
	double check_cost;

	do
	{
		check_cost = (*best_cost);

		int i = 0;
		for (; i < tsp_in->num_nodes; i++)//aggiunto i=0 dentro
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
					double cost_j_h=0.0; //cost[i, succ[i]]
					double cost_i_j=0.0;
					double cost_k_h=0.0;
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
						//dist(i, succ[i], tsp_in, &cost_i_k);//aggiunto
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
	} while (check_cost != (*best_cost)); //abs(check_cost - (*cost))<1e-10

	visited_nodes[0] = 0;
	int i = 0;
	int count = 1;
	int node = -1;
	for (; count < tsp_in->num_nodes; count++)
	{
		node = succ[i];
		i = node;
		visited_nodes[count] = node;  //inserirlo nell'algoritmo
		//printf("[ %d ] count: %d\n", greedy_count, visited_nodes[count]);
	}

	free(succ);
}

void vns(tsp_instance* tsp_in, int* visited_nodes, double *best_cost, double deadline)
{
	srand(time(NULL));

	double remaining_time = deadline;

	int* local_min_visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
		local_min_visited_nodes[i] = visited_nodes[i];

	int num_local_mins = 0;
	int num_iterations = 0;
	double local_min_cost = (*best_cost);


	while(remaining_time > 0 && num_local_mins < MAX_LOCAL_MINS && num_iterations < MAX_NUM_ITERATIONS)
	{
		time_t start = clock();

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

		time_t end = clock();
		remaining_time = remaining_time - ((double)(end - start) / CLOCKS_PER_SEC);

		#ifndef MULTI_START
			printf("%sRemaining time : %s%.2lf\n", CYAN, WHITE, remaining_time);
		#endif 

	}
}

void update_solution(int* visited_nodes, double* sol, int num_nodes)//necessario aver già allocato sol , un vettore di double generico
{
	int num_edges = num_nodes * (num_nodes - 1) / 2;

	int i;
	for (i = 0; i < num_edges; i++)
		sol[i] = 0.0;

	for (i = 0; i < num_nodes; i++)
		sol[generic_xpos(visited_nodes[i], visited_nodes[(i + 1) % num_nodes], num_nodes)] = 1.0;

}

void succ_construction(int* visited_nodes, int* succ, int num_nodes)//succ deve essere gia' allocato
{
	int i;
	for (i = 0; i < num_nodes; i++)
		succ[visited_nodes[i]] = visited_nodes[(i+1) % num_nodes];

}

void tabu_search(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, double deadline)
{
	time_t start = clock();
	double remaining_time = deadline;

	int* succ = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	succ_construction(visited_nodes, succ, tsp_in->num_nodes);

	#ifndef MULTI_START
	if (tsp_in->verbose > 50)
		printf("%sStarting cost :%s%.2lf\n",RED, WHITE, *best_cost);
	#endif 

	/*greedy_refinement(tsp_in, visited_nodes, best_cost); //prima di iniziare l'algoritmo devo essere sicura di trovarmi in un minimo locale

	#ifndef MULTI_START
	if (tsp_in->verbose > 50)
		printf("%sCost after first greedy :%s %.2lf\n",GREEN, WHITE, *best_cost);
	#endif 
	*/
	int min_tenure = ceil(tsp_in->num_nodes / 10.0);
	int max_tenure = ceil(tsp_in->num_nodes / 5.0);

	int** tabu_list = (int**)calloc((size_t)2, sizeof(int*));//[tabu_list[0][i], tabu_list[1][i]]
	tabu_list[0] = (int*)calloc((size_t)max_tenure, sizeof(int));
	tabu_list[1] = (int*)calloc((size_t)max_tenure, sizeof(int));

	tabu_list_params param;
	param.end_list = -1;
	param.start_list = 0;

	int i;
	for (i = 0; i < max_tenure; i++)
	{
		tabu_list[0][i] = -1;
		tabu_list[1][i] = -1;
	}

	int num_tabu_edges = 0;
	double actual_cost = *best_cost;

	time_t end = clock();
	remaining_time = remaining_time - ((double)(end - start) / CLOCKS_PER_SEC);

	int num_iteration = 0;
	for (; num_iteration < MAX_NUM_ITERATIONS && remaining_time > 0; num_iteration++)
	{
		start = clock();

		double min_increase = DBL_MAX;

		min_increase = move2opt_for_tabu_search(tsp_in, succ, tabu_list, &max_tenure , &param);

		if (num_tabu_edges < max_tenure)
			num_tabu_edges = num_tabu_edges + 2;

		actual_cost = actual_cost + min_increase;

		if (min_increase < 0.0 )
		{
			greedy_refinement_for_tabu_search(tsp_in, succ, tabu_list, &param, max_tenure, min_tenure, &num_tabu_edges, &actual_cost);

			#ifndef MULTI_START
			if (tsp_in->verbose > 50)
				printf("%sfind local minimum with cost : %s%.2lf\n",GREEN, WHITE, actual_cost);
			#endif 

		}//terminato il greedy devo riampliare la tabu list fino al massimo
		else
		{
			#ifndef MULTI_START
			if (tsp_in->verbose > 50)
				printf("%sUpdate cost :%s %.2lf\n", BLUE, WHITE, actual_cost);
			#endif
		}
	
		if (actual_cost < *best_cost)//controllare precisione
		{
			*best_cost = actual_cost;
			
			int j = 0;
			int next = j;
			for (; j < tsp_in->num_nodes; j++)
			{
				visited_nodes[j] = next;
				next = succ[next];
			}
		}

		end = clock();
		remaining_time = remaining_time - ((double)(end - start) / CLOCKS_PER_SEC);

#ifndef MULTI_START
		printf("%sRemaining time : %s%.2lf\n", CYAN, WHITE, remaining_time);
#endif 

	}


	free(succ);
	free(tabu_list[0]);
	free(tabu_list[1]);
	free(tabu_list);
	
}

double move2opt_for_tabu_search(tsp_instance* tsp_in, int* succ, int** tabu_list, int* tenure, tabu_list_params* params) //return the delta_min
{
	double delta_min = DBL_MAX;
	int start_edge1 = -1;
	int end_edge1 = -1;
	int start_edge2 = -1;
	int end_edge2 = -1;

	int i=0;
	for (; i < tsp_in->num_nodes; i++)//aggiunto i=0 dentro
	{
		double cost_i_k = 0.0; //cost [i, succ[i]]

		int j = 0;
		for (; j < tsp_in->num_nodes; j++)
		{

			if (j != i && j != succ[i] && succ[j] != i && succ[j] != succ[i])
			{

				if (check_tabu_list(tabu_list, *tenure, i, j) == 1)
					continue;

				if (check_tabu_list(tabu_list, *tenure, succ[i], succ[j]) == 1)
					continue;


				double cost_j_h = 0.0; //cost[i, succ[i]]
				double cost_i_j = 0.0;
				double cost_k_h = 0.0;
				if (tsp_in->integerDist)
				{
					int x1, x2, x3, x4;
					dist(j, succ[j], tsp_in, &x1);
					dist(i, j, tsp_in, &x2);
					dist(succ[i], succ[j], tsp_in, &x3);
					dist(i, succ[i], tsp_in, &x4);
					cost_j_h = (double)x1;
					cost_i_j = (double)x2;
					cost_k_h = (double)x3;
					cost_i_k = (double)x4;
				}
				else
				{
					dist(j, succ[j], tsp_in, &cost_j_h);
					dist(i, j, tsp_in, &cost_i_j);
					dist(succ[i], succ[j], tsp_in, &cost_k_h);
					dist(i, succ[i], tsp_in, &cost_i_k);//aggiunto
				}
				double delta = cost_i_j + cost_k_h - cost_i_k - cost_j_h;

				if ( fabs(delta) > 0.1   && delta < delta_min)
				{
					delta_min = delta;
					start_edge1 = i;
					end_edge1 = j;
					start_edge2 = succ[i];
					end_edge2 = succ[j];
				}
			}
		}
	}

	add_element(tabu_list[0], tabu_list[1], *tenure, start_edge1, start_edge2, 0, 0, params);
	add_element(tabu_list[0], tabu_list[1], *tenure, end_edge1, end_edge2, 0, 0, params);

	int k = start_edge2;
	int next = start_edge2;

	int count = 0;
	int* orientation = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	while (next != end_edge1)
	{
		orientation[count++] = next;
		next = succ[next];
	}
	next = succ[k];
	int next_count = 0;

	while (next != end_edge1)
	{
		int temp = succ[next];
		succ[next] = orientation[next_count++];
		next = temp;
	}
	succ[start_edge1] = end_edge1;
	succ[k] = succ[end_edge1];
	succ[end_edge1] = orientation[count - 1];
	free(orientation);

	return delta_min;
}

int check_tabu_list(int** tabu_list, int tenure, int node1, int node2)//return 1 if the edge is forbidden, 0 otherwise
{
	int i;
	for (i = 0; i < tenure; i++)
	{
		if ((tabu_list[0][i] == node1 && tabu_list[1][i] == node2) ||
			(tabu_list[0][i] == node2 && tabu_list[1][i] == node1)    )
		{
			return 1;
		}

	}

	return 0;
}

void add_element(int* list1, int* list2, int dimension, int element1, int element2, int with_reduction, int logically_full, tabu_list_params* param)
{
	if (!with_reduction)
	{
		param->end_list = (param->end_list + 1) % dimension;
		list1[param->end_list] = element1;
		list2[param->end_list] = element2;

		if (param->end_list == param->start_list || logically_full)//the buffer il full
			param->start_list = (param->start_list + 1) % dimension;
	}
	else
	{
		list1[param->start_list] = -1;
		list2[param->start_list] = -1;

		list1[(param->start_list + 1) % dimension] = -1;
		list2[(param->start_list + 1) % dimension] = -1;

		param->start_list = (param->start_list + 2) % dimension;

		param->end_list = (param->end_list + 1) % dimension;

		list1[param->end_list] = element1;
		list2[param->end_list] = element2;
	}
}

void greedy_refinement_for_tabu_search(tsp_instance* tsp_in, int* succ, int** tabu_list, tabu_list_params* param, int max_tenure, 
	int min_tenure, int* num_tabu_edges, double* best_cost)
{
	//la tabu list ad ogni iterazione deve diminuire la sua dimensione, prima la scelta e poi la riduzione
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
							#ifdef REACTIVE
								if (*num_tabu_edges > min_tenure)
								{
									add_element(tabu_list[0], tabu_list[1], max_tenure, i, succ[i], 1, 0, param);
									(*num_tabu_edges)--;
									add_element(tabu_list[0], tabu_list[1], max_tenure, j, succ[j], 1, 0, param);
									(*num_tabu_edges)--;
								}
								else if (*num_tabu_edges == min_tenure)
								{
									add_element(tabu_list[0], tabu_list[1], max_tenure, i, succ[i], 0, 1, param);
									add_element(tabu_list[0], tabu_list[1], max_tenure, j, succ[j], 0, 1, param);
								}
								else
								{
									add_element(tabu_list[0], tabu_list[1], max_tenure, i, succ[i], 0, 0, param);
									add_element(tabu_list[0], tabu_list[1], max_tenure, j, succ[j], 0, 0, param);
								
									*num_tabu_edges = *num_tabu_edges + 2;
								}
						
							#else

								add_element(tabu_list[0], tabu_list[1], max_tenure, i, succ[i], 0, 0, param);
								add_element(tabu_list[0], tabu_list[1], max_tenure, j, succ[j], 0, 0, param);	

								if (*num_tabu_edges < max_tenure)
									*num_tabu_edges = *num_tabu_edges + 2;

							#endif

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
}


void simulated_annealing(tsp_instance* tsp_in, int* visited_nodes, double* best_cost, double deadline)
{
	//visited_nodes contiene la prima soluzione calcolata esternamente alla funzione, best cost il suo costo
	#ifndef MULTI_START

	printf("%sStarting cost:%s %.2lf\n", RED, WHITE, *best_cost);

	#endif 

	int outer_iteration;

	double remaining_time = deadline;

	double alpha = 0.99;
	
	double t_min = 100.0; 

	double t_max = 5000.0;

	double t = t_max ;

	double x = 2.3; // 10 = e^x

	double r = 0.0;

	int* new_visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

	int k;
	for (k = 0; k < tsp_in->num_nodes; k++)
		new_visited_nodes[k] = visited_nodes[k];

	double cost = *best_cost;
	
	//formula for the decrease of the temperature: t = (alpha)^i * t_start + t_min

	for (outer_iteration = 0; t - t_min > 0.1 && remaining_time > 0.0 ; outer_iteration++)
	{
		int i;
		int increase_accepted = 0;

#ifndef MULTI_START

		printf("%stemperature :%s %.2lf   ", BLUE, WHITE, t);

#endif 

		for (i=0; !increase_accepted; i++)
		{
			srand( (outer_iteration+1) * (i+1) * 100);

			time_t start_inner_iteration = clock();

			double choice =(double) (rand() % 100);
			double range = 100.0 / tsp_in->num_nodes;
			int j;
			for (j = 0; j < tsp_in->num_nodes; j++)
			{
				if (choice < ((j+1) * range) )
					break;
			}

			int index_node1 = -1;
			if (j != tsp_in->num_nodes)
				index_node1 = j;
			else
				printf("ERROR\n");

			int value = rand() % 100;
			int offset = (value < 50) ?	1 : 2;

			int  index_node2 = (index_node1 + offset) %tsp_in->num_nodes;

			//printf("\033[0;34m index_node1 = %d     index_node2=%d\n", index_node1, index_node2);
			//printf("\033[0m");

			double new_cost = 0.0;
		
			double delta = 0.0;
			if (tsp_in->integerDist)
			{
				int  x1_new, x2_new, x1_old, x2_old;
				dist(new_visited_nodes[(index_node1 - 1 + tsp_in->num_nodes)% tsp_in->num_nodes], new_visited_nodes[index_node2], tsp_in, &x1_new);
				dist(new_visited_nodes[index_node1], new_visited_nodes[(index_node2  + 1) % tsp_in->num_nodes], tsp_in, &x2_new);
				dist(new_visited_nodes[(index_node1 - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], new_visited_nodes[index_node1], tsp_in, &x1_old);
				dist(new_visited_nodes[index_node2], new_visited_nodes[(index_node2 + 1)% tsp_in->num_nodes], tsp_in, &x2_old);
				delta = (double)(x1_new + x2_new - x1_old - x2_old);
			}
			else
			{
				double  x1_new, x2_new, x1_old, x2_old;
				dist(new_visited_nodes[(index_node1 - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], new_visited_nodes[index_node2], tsp_in, &x1_new);
				dist(new_visited_nodes[index_node1], new_visited_nodes[(index_node2 + 1) % tsp_in->num_nodes], tsp_in, &x2_new);
				dist(new_visited_nodes[(index_node1 - 1 + tsp_in->num_nodes) % tsp_in->num_nodes], new_visited_nodes[index_node1], tsp_in, &x1_old);
				dist(new_visited_nodes[index_node2], new_visited_nodes[(index_node2 + 1) % tsp_in->num_nodes], tsp_in, &x2_old);
				delta = x1_new + x2_new - x1_old - x2_old;
			}

			new_cost = cost + delta;

			//printf("new cost : %lf <<<<<<<<<<< delta %lf >>>>>>>>> best cost : %lf\n", new_cost,delta, *best_cost);

			if (new_cost < cost)
			{
				/*
				printf("\033[1;32m Improvement  cost: %lf\n", new_cost);
				printf("node1: visited_nodes[%d]=%d     node2: visited_nodes[%d] = %d\n",
					index_node1, new_visited_nodes[index_node1], index_node2, new_visited_nodes[index_node2]);
				printf("\033[0m");
				*/

				#ifndef MULTI_START

					printf("%simprovement cost:%s %.2lf\n", GREEN, WHITE, new_cost);

				#endif 

				int temp = new_visited_nodes[index_node1];
				new_visited_nodes[index_node1] = new_visited_nodes[index_node2];
				new_visited_nodes[index_node2] = temp;

				cost = new_cost;

				greedy_refinement(tsp_in, new_visited_nodes, &cost);

				if (cost < *best_cost)
				{
					int h;
					for (h = 0; h < tsp_in->num_nodes; h++)
						visited_nodes[h] = new_visited_nodes[h];

					*best_cost = cost;
				}
				
			}
			else if (new_cost > cost)
			{
				int accept = 0;
				
				double expo = (new_cost - cost) / t; //esponende della e nel calcolo della probabilità
				int magnitude =(int) floor(expo / x); //ordine di grandezza della probabilità
				double coeff = exp(expo - magnitude * x);

				int k;
				int count = 0;
				for (k = 0; k < magnitude; k++)
				{
					int r = rand() % 100;
					if (r < 10)
						count++;
				}
				if (count == magnitude)
				{
					int max_range = (int) (10 * coeff);
					
					int r = rand() % max_range;
					if (r < 10)
						accept = 1;
				}
				//printf("\033[1;36m magnitude = %d    coeff = %lf \n", magnitude, coeff);
				//printf("\033[0m");
				
				if (accept)
				{
					/*
					printf("\033[1;31m Accepted  new solution       cost: %lf\n",new_cost);
					printf("node1: visited_nodes[%d]=%d     node2: visited_nodes[%d] = %d\n",
						index_node1, new_visited_nodes[index_node1], index_node2, new_visited_nodes[index_node2]);
					printf("\033[0m");
					*/

					#ifndef MULTI_START

						printf("%snew worst cost:%s   %.2lf\n", BLUE, WHITE, new_cost);

					#endif 

					

					int temp = new_visited_nodes[index_node1];
					new_visited_nodes[index_node1] = new_visited_nodes[index_node2];
					new_visited_nodes[index_node2] = temp;

					cost = new_cost;

					increase_accepted = 1;
				}

			}

			time_t end_inner_iteration = clock();

			remaining_time = remaining_time - ((double)(end_inner_iteration - start_inner_iteration) / CLOCKS_PER_SEC);

			if (remaining_time < 0)
				break;
		}
		t = ((pow(alpha, outer_iteration+1) )* t_max) + t_min;
	}
}
void genetic_solver(tsp_instance* tsp_in)
{
	time_t start = clock();

	tsp_in->bestCostD = DBL_MAX;
	tsp_in->bestCostI = INT_MAX;

	int num_worst = 0;
	int* worst_members = (int*)calloc((size_t)NUM_WORST_MEMBERS, sizeof(int));

	//printf(LINE);
	printf("%sHeuristic solver%s\n", RED, WHITE);
	printf("%s[Construction]  %sNearest Neighborhood\n", BLUE, WHITE);

	printf("%s[Meta-heuristic] %sGenetic\n", BLUE, WHITE);
	printf(LINE);

	int** members = (int**)calloc((size_t)POPULATION_SIZE, sizeof(int*));
	double* fitnesses = (double*)calloc((size_t)POPULATION_SIZE, sizeof(double));
	double sum_prob = 0.0;

	int i = 0;
	for (i = 0; i < POPULATION_SIZE; i++)
	{
		members[i] = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		fitnesses[i] = DBL_MAX;
	}

	int num_threads = NUM_MULTI_START;
	pthread_t* threads = (pthread_t*)calloc((size_t)num_threads, sizeof(pthread_t));
	construction_args* param = (construction_args*)calloc((size_t)num_threads, sizeof(construction_args));

	int update = 0;

	printf("[");

	i = 0;
	for (; i < 10; i++)
		printf(" ");

	printf("] %3d %%     ", update);

	if (tsp_in->integerDist)
	{
		int worst = 0;
		printf("worst cost: INF     ");
		printf("incumbent: INF     ");
		printf("average: INF");
	}
	else
	{
		double worst = 0.0;
		printf("worst cost: INF     ");
		printf("incumbent: INF     ");
		printf("average: INF");
	}

	int num_instances = 0;
	int num_members = POPULATION_SIZE / num_threads;
	int best_index = 0;
	double sum_fitnesses = 0.0;

	i = 0;
	for (; i < num_threads; i++)
	{
		param[i].tsp_in = tsp_in;
		param[i].members = members;
		param[i].fitnesses = fitnesses;
		param[i].num_instances = &num_instances;
		param[i].first_index = (POPULATION_SIZE / num_threads) * i;
		param[i].best_index = &best_index;
		param[i].sum_fitnesses = &sum_fitnesses;

		if (i == num_threads)
			param[i].num_members = num_members;
		else
		{
			num_members += POPULATION_SIZE % num_threads;
			param[i].num_members = num_members;
		}

		param[i].sum_prob = &sum_prob;
		pthread_create(&(threads[i]), NULL, construction, (void*)&param[i]);
	}

	for (i = 0; i < num_threads; i++)
	{
		int rc = pthread_join(threads[i], NULL);

		if (rc)
			exit(-1);
	}

	free(threads);
	free(param);

	for (i = 0; i < NUM_WORST_MEMBERS; i++)
	{
		worst_members[i] = -1;
	}
	/*
	int* succ = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int n_comps = 1;*/
	int  num_edges = (tsp_in->num_nodes * (tsp_in->num_nodes-1)) / 2;
	tsp_in->sol = (double*)calloc((size_t)num_edges, sizeof(double));

	/*update_solution(members[best_index], tsp_in->sol, tsp_in->num_nodes);
	define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
	plot(tsp_in, succ, comp, &n_comps);

	double tmp = 0.0;
	for (i = 0; i < tsp_in->num_nodes; i++)
	{
		double x;
		dist(members[best_index][i], members[best_index][(i + 1) % tsp_in->num_nodes], tsp_in, &x);
		tmp += x;
	}*/
	evolution(tsp_in, members, fitnesses, &best_index, worst_members, &sum_prob, &sum_fitnesses, start);

	print_cost(tsp_in);
	printf("%sExecution time:%s %.2lf seconds\n", GREEN, WHITE, tsp_in->execution_time);
	printf("%s%s%s", RED, LINE, WHITE);

	
	int* succ = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int* comp = (int*)calloc(tsp_in->num_nodes, sizeof(int));
	int n_comps = 1;
	
	succ_construction(members[best_index], succ, tsp_in->num_nodes);
	update_solution(members[best_index], tsp_in->sol, tsp_in->num_nodes);

	if (tsp_in->plot)
	{
		int k;
		for (k = 0; k < tsp_in->num_nodes; k++)
			comp[k] = 1;

		plot(tsp_in, succ, comp, &n_comps);

		free(succ);
		free(comp);
	}

	for (i = 0; i < POPULATION_SIZE; i++)
		free(members[i]);

	free(members);
	free(fitnesses);
	free(tsp_in->sol);
}

void construction(void* param)
{
	construction_args* args = (construction_args*)param;
	double best_cost = DBL_MAX;
	int best_index = -1;
	double sum_prob = 0.0;
	double sum_fitnesses = 0.0;
	double best_wrong_cost=0.0;

	int i = 0;
	for (; i < args->num_members && i < POPULATION_SIZE; i++)
	{
		nearest_neighborhood(args->tsp_in, args->members[i + args->first_index], &(args->fitnesses[i + args->first_index]), args->first_index + i + 1, (i + args->first_index) % args->tsp_in->num_nodes);
		double cost = args->fitnesses[i+args->first_index];
		greedy_refinement(args->tsp_in, args->members[i + args->first_index], &(args->fitnesses[i + args->first_index]));
		
		sum_prob += (1000.0 / args->fitnesses[i + args->first_index]);
		sum_fitnesses += args->fitnesses[i + args->first_index];

		if ((args->fitnesses[i + args->first_index]) < best_cost)
		{
			best_cost = args->fitnesses[i + args->first_index];
			best_index = i + args->first_index;
			best_wrong_cost = cost;
		}
	}

	pthread_mutex_lock(&mutex);

	*(args->num_instances) += args->num_members;
	/*printf("\r[%s", GREEN);

	i = 0;
	int update = (((double)*(args->num_instances)) / (double)(POPULATION_SIZE)) * 10.0;

	for (; i < 10; i++)
	{
		if ((double)i < update)
			printf("=");
		else
			printf(" ");
	}

	printf("%s] %3d %%     ", WHITE, (int)(update * 10.0));
	*/

	if (args->tsp_in->integerDist)
	{
		if (best_cost < args->tsp_in->bestCostI)
		{
			args->tsp_in->bestCostI = (int)best_cost;
			*args->best_index = best_index;
		}

		/*
		printf("incumbent: %d     ", args->tsp_in->bestCostI);
		printf("best_index: %d", (*args->best_index));
		printf("average: %d", (int)(*args->sum_fitnesses / (double)POPULATION_SIZE));
		*/
		/*
		printf("[Thread %2d]  best_cost: %d    instances: %5d     incumbent: %d\n",
			   args->id, best_cost, *(args->num_instances), args->tsp_in->bestCostI);
		*/
	}
	else
	{
		if (best_cost < args->tsp_in->bestCostD)
		{
			args->tsp_in->bestCostD = best_cost;
			*args->best_index = best_index;
			printf("[%d] COMPARE %.2lf   wrong: %.2lf", best_index, best_cost, best_wrong_cost);
		}

		/*
		printf("incumbent: %.2lf     ", args->tsp_in->bestCostD);
		printf("best_index: %d     ", (*args->best_index));
		printf("average: %.2lf", *args->sum_fitnesses / (double)POPULATION_SIZE);
		*/

		/*
		printf("[Thread %2d]  best_cost: %.2lf    instances: %5d     incumbent: %.2lf\n",
			   args->id, best_cost, *(args->num_instances), args->tsp_in->bestCostD);
		*/
	}

	(*(args->sum_prob)) += sum_prob;
	(*(args->sum_fitnesses)) += sum_fitnesses;

	pthread_mutex_unlock(&mutex);

	pthread_exit(NULL);
}

void evolution(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index, 
	           int* worst_members, double* sum_prob, double* sum_fitnesses, time_t start)
{
	int num_epochs = 0;
	time_t end = clock();
	printf("\n%s", LINE);
	printf("%s[Construction phase]  ", RED);
	printf("%sbest index:%s %5d\n", GREEN, WHITE, *best_index);
	int index=0; //first_index of worst_members

	while (num_epochs < 150 && (((double)(end - start) / (double)CLOCKS_PER_SEC)) < tsp_in->deadline)
	{
		//printf("Sum_fitnesses: %.2lf      \n", (*sum_fitnesses));
		if (index%NUM_WORST_MEMBERS == 0)
		{
			index = 0;
			update_worst(fitnesses, worst_members);
		}

		if (num_epochs % 5 == 0)
		{
			crossover(tsp_in, members, fitnesses, best_index, worst_members, sum_prob, sum_fitnesses, (num_epochs + 1) * 100, &index);

			printf(LINE);
			if (tsp_in->integerDist)
				printf("%s[Crossover]%s      added instances: %d     incumbent: %d    average: %d\n",
					RED, WHITE, index, tsp_in->bestCostI, (int)(*sum_fitnesses / (double)POPULATION_SIZE));
			else
				printf("%s[Crossover]%s      added instances: %d     incumbent: %.2lf    average: %.2lf\n",
					RED, WHITE, index, tsp_in->bestCostD, *sum_fitnesses / (double)POPULATION_SIZE);
		}
		else
		{
			mutation(tsp_in, members, fitnesses, best_index, worst_members, sum_prob, sum_fitnesses, (num_epochs + 1) * 100, &index);

			if (tsp_in->integerDist)
				printf("%s[Mutation]%s       added instances: %d     incumbent: %d    average: %d\n", 
					   BLUE, WHITE, index, tsp_in->bestCostI, (int) (*sum_fitnesses/(double) POPULATION_SIZE));
			else
				printf("%s[Mutation]%s       added instances: %d     incumbent: %.2lf    average: %.2lf\n",
					BLUE, WHITE, index, tsp_in->bestCostD, *sum_fitnesses / (double)POPULATION_SIZE);
		}

		printf("%sbest index:%s %5d\n", GREEN, WHITE, *best_index);

		end = clock();
		num_epochs++;
	}

	printf(LINE);
	tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
}

void crossover(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index,
	          int* worst_members, double* sum_prob, double* sum_fitnesses, int seed, int* index)
{
	srand(seed);

	int i = 0;
	for (; i < NUM_WORST_MEMBERS/2; i++)
	{
		int j = 0;
		double dad = (double)(rand() % ((int)((*sum_prob) * 100000.0)));
		double mom = (double)(rand() % ((int)((*sum_prob) * 100000.0)));
		//printf("SUM: %.2lf\n",sum_prob*1000);
		//printf("SUM: %.2lf\n", dad);
		//printf("SUM: %.2lf\n", mom);

		double sum_ranges = 0.0;
		int found_dad = 0;
		int found_mom = 0;
		int dad_index;
		int mom_index;

		/*
		printf(LINE);
		printf("Mom: %.2lf   fitness: %.2lf\n", mom, *sum_fitnesses/ tsp_in->num_nodes);
		printf("Dad: %.2lf   fitness: %.2lf\n", dad, *sum_fitnesses/ tsp_in->num_nodes);
		*/

		//1/fitness[i]
		for (j = 0; j < POPULATION_SIZE; j++)
		{

			if (((100000000.0 / ((*sum_prob) * fitnesses[j])) + sum_ranges) > dad && !found_dad)
			{
				dad_index = j;
				found_dad = 1;
			}

			if (((100000000.0 / ((*sum_prob) * fitnesses[j])) + sum_ranges) > mom && !found_mom)
			{
				mom_index = j;
				found_mom = 1;
			}

			if (found_dad && found_mom)
				break;

			sum_ranges += (100000000.0 / ((*sum_prob) * fitnesses[j]));
		}

		int* offspring1 = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		int* offspring2 = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));

		j = 0;
		int begin = (int)(((double)tsp_in->num_nodes) /2.0);

		for (j = begin; j < tsp_in->num_nodes; j++)
		{
			offspring1[j] = members[dad_index][j];
			offspring2[j-begin] = members[mom_index][j-begin];
		}

		int count1 = 0;
		int count2 = tsp_in->num_nodes - begin;
		for (j = 0; j < tsp_in->num_nodes && (count1 < begin || count2 < tsp_in->num_nodes); j++)
		{
			int k;

			//Offspring 1 fullfill
			for (k = begin; k < tsp_in->num_nodes; k++)
			{
				if (members[mom_index][j] == offspring1[k])
					break;
			}

			if (k == tsp_in->num_nodes)
				offspring1[count1++] = members[mom_index][j];


			//Offspring2 fullfill
			for (k = 0; k < (tsp_in->num_nodes - begin); k++)
			{
				if (members[dad_index][j] == offspring2[k])
					break;
			}

			if (k == (tsp_in->num_nodes-begin))
				offspring2[count2++] = members[dad_index][j];
		}

		int k;
		double fitness = 0.0;
		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			members[worst_members[(*index)]][k] = offspring1[k];

			if (tsp_in->integerDist)
			{
				int x;
				dist(offspring1[k], offspring1[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += (double)x;
			}
			else
			{
				double x;
				dist(offspring1[k], offspring1[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += x;
			}
		}

		free(offspring1);
		greedy_refinement(tsp_in, members[worst_members[*index]], &fitness);

		if (tsp_in->integerDist)
		{
			if ((int)fitness < tsp_in->bestCostI)
			{
				tsp_in->bestCostI = (int)fitness;
				(*best_index) = worst_members[(*index)];
			}
		}
		else
		{
			if (fitness < tsp_in->bestCostD)
			{
				tsp_in->bestCostD = fitness;
				(*best_index) = worst_members[(*index)];
			}
		}

		(*sum_prob) += ((1000.0 / fitness) - (1000.0 / fitnesses[worst_members[(*index)]]));
		(*sum_fitnesses) += (fitness - fitnesses[worst_members[(*index)]]);
		fitnesses[worst_members[(*index)]] = fitness;
		worst_members[(*index)] = -1;
		(*index)++;

		fitness = 0.0;
		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			members[worst_members[(*index)]][k] = offspring2[k];

			if (tsp_in->integerDist)
			{
				int x;
				dist(offspring2[k], offspring2[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += (double)x;
			}
			else
			{
				double x;
				dist(offspring2[k], offspring2[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += x;
			}
		}

		free(offspring2);
		greedy_refinement(tsp_in, members[worst_members[*index]], &fitness);

		if (tsp_in->integerDist)
		{
			if ((int)fitness < tsp_in->bestCostI)
			{
				tsp_in->bestCostI = (int)fitness;
				(*best_index) = worst_members[(*index)];
			}
		}
		else
		{
			if (fitness < tsp_in->bestCostD)
			{
				tsp_in->bestCostD = fitness;
				(*best_index) = worst_members[(*index)];
			}
		}

		(*sum_prob) += ((1000.0 / fitness) - (1000.0 / fitnesses[worst_members[(*index)]]));
		(*sum_fitnesses) += (fitness - fitnesses[worst_members[(*index)]]);
		fitnesses[worst_members[(*index)]] = fitness;
		worst_members[(*index)] = -1;
		(*index)++;

	}
	//printf("Index: %d   ", *index);
}

void mutation(tsp_instance* tsp_in, int** members, double* fitnesses, int* best_index,
	          int* worst_members, double* sum_prob, double* sum_fitnesses, int seed, int* index)
{
	srand(seed);

	int i = 0;
	for (; i < NUM_WORST_MEMBERS; i++)
	{
		int j = 0;
		double dad = (double)(rand() % ((int)((*sum_prob) * 100000.0)));
		int dad_index;

		double sum_ranges = 0.0;

		//1/fitness[i]
		for (j = 0; j < POPULATION_SIZE; j++)
		{
			if (((100000000.0 / ((*sum_prob) * fitnesses[j])) + sum_ranges) > dad)
			{
				dad_index = j;
				break;
			}

			sum_ranges += (100000000.0 / ((*sum_prob) * fitnesses[j]));
		}

		double start_choice = (double)(rand() % 100);
		double end_choice = (double)(rand() % 100);
		int start_range = -1;
		int end_range = -1;
		double node_range = 200.0 / tsp_in->num_nodes;
		int found_start = 0;
		int found_end = 0;

		for (j = 0; j < tsp_in->num_nodes; j++)
		{
			if ((node_range * (j + 1)) > start_choice && !found_start)
			{
				start_range = j;
				found_start = 1;
			}

			if ((node_range * (j + 1)) > end_choice && !found_end)
			{
				end_range = j + (int)((double)tsp_in->num_nodes / 2);
				found_end = 1;
			}

			if (found_start && found_end)
				break;

		}

		int* offspring = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
		j = 0;
		start_range = (int)(((double)tsp_in->num_nodes) / 2.0);
		end_range = (int)(((double)tsp_in->num_nodes));

		for (j = start_range; j < end_range; j++)
		{
			offspring[j] = members[dad_index][end_range - 1 - j + start_range];
		}

		for (j = 0; j < start_range; j++)
		{
			offspring[j] = members[dad_index][j];
		}

		for (j = end_range; j < tsp_in->num_nodes; j++)
		{
			offspring[j] = members[dad_index][j];
		}

		int k;
		double fitness = 0.0;

		for (k = 0; k < tsp_in->num_nodes; k++)
		{
			members[worst_members[(*index)]][k] = offspring[k];

			if (tsp_in->integerDist)
			{
				int x;
				dist(offspring[k], offspring[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += (double)x;
			}
			else
			{
				double x;
				dist(offspring[k], offspring[(k + 1) % tsp_in->num_nodes], tsp_in, &x);

				fitness += x;
			}
		}

		free(offspring);
		greedy_refinement(tsp_in, members[worst_members[*index]], &fitness);

		if (tsp_in->integerDist)
		{
			if ((int)fitness < tsp_in->bestCostI)
			{
				tsp_in->bestCostI = (int)fitness;
				(*best_index) = worst_members[(*index)];
			}
		}
		else
		{
			if (fitness < tsp_in->bestCostD)
			{
				tsp_in->bestCostD = fitness;
				(*best_index) = worst_members[(*index)];
			}
		}

		(*sum_prob) += ((1000.0 / fitness) - (1000.0 / fitnesses[worst_members[(*index)]]));
		(*sum_fitnesses) += (fitness - fitnesses[worst_members[(*index)]]);
		fitnesses[worst_members[(*index)]] = fitness;
		worst_members[(*index)] = -1;
		(*index)++;
	}
	//printf("Index: %d   ", *index);
}

void update_worst(double* fitnesses, int* worst_members)
{
	int i = 0;
	int size = 0;

	for (; i < POPULATION_SIZE; i++)
	{
		int j = 0;
		int found_worst = 0;

		for (; j <= size; j++)
		{
			if (worst_members[j] == -1)
				break;

			if (fitnesses[worst_members[j]] < fitnesses[i])
			{
				found_worst = 1;
				break;
			}
		}

		if (j < NUM_WORST_MEMBERS)
		{
			int k = j;

			for (j = size; found_worst && j > k; j--)
				worst_members[j] = worst_members[j - 1];

			worst_members[j] = i;
			size = (size == (NUM_WORST_MEMBERS - 1)) ? size : (size + 1);
		}
	}

	/*
	for(i=0; i<1000; i++)
		printf("[%d] worst_member: %2d    fitness: %.2lf\n", i, worst_members[i], fitnesses[worst_members[i]]);
	*/
}