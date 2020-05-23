#include "heuristic.h"
#include "utility.h"

void heuristic_solver(tsp_instance* tsp_in)
{
	int* visited_nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	int num_edges = tsp_in->num_nodes * (tsp_in->num_nodes - 1) / 2;
	tsp_in->sol =(double *) calloc( num_edges ,sizeof(double));

	switch (tsp_in->alg)
	{
		case 6:
		{
			printf(LINE);
			printf("Heuristic nearest neighborhood solver\n");
			time_t start = clock();
			nearest_neighborhood(tsp_in, visited_nodes);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}
		case 7:
		{
			printf(LINE);
			printf("Heuristic insertion solver\n");
			time_t start = clock();
			insertion(tsp_in, visited_nodes);
			time_t end = clock();
			tsp_in->execution_time = ((double)(end - start) / (double)CLOCKS_PER_SEC);
			break;
		}
	}
	
	update_solution(visited_nodes, tsp_in->sol, tsp_in->num_nodes);
	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	int* succ =(int *) calloc(tsp_in->num_nodes, sizeof(int));
	int* comp =(int *) calloc(tsp_in->num_nodes, sizeof(int));
	int n_comps = 1;
	if (tsp_in->plot)
	{
		switch (tsp_in->alg)
		{
		case 6: case 7:
		{
			//define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
			succ_construction(visited_nodes, succ, tsp_in->num_nodes);
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

	//define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);
	double cost;
	greedy_refinement(tsp_in, visited_nodes, &cost);
	update_solution(visited_nodes, tsp_in->sol, tsp_in->num_nodes);//solo se soluzione migliorata

	print_cost(tsp_in);
	printf("Execution time: %lf\n", tsp_in->execution_time);

	succ = calloc(tsp_in->num_nodes, sizeof(int));
	comp = calloc(tsp_in->num_nodes, sizeof(int));
	n_comps = 1;
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

void nearest_neighborhood(tsp_instance* tsp_in, int* visited_nodes)
{
	int i = 0;

	int* nodes = (int*)calloc((size_t)tsp_in->num_nodes, sizeof(int));
	//vettore di dim num nodes -> nodes utilizzato per fare le verifiche in min_cost
	//metto a 1 i nodi visitati

	if (tsp_in->integerDist)
		tsp_in->bestCostI = 0;
	else
		tsp_in->bestCostD = 0.0;

	visited_nodes[0] = 0;
	nodes[0] = 1;
	int count = 1;

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
				min_dist = (double) x;
			}
			else
				dist(i, 0, tsp_in, &min_dist);

			printf("count = %d\n", count);
			
		}
		else
		{
			visited_nodes[count] = best;//aggiungere nodo a visited nodes
			count++;
		}

		if (tsp_in->integerDist)
			tsp_in->bestCostI += (int) min_dist;
		else
			tsp_in->bestCostD += min_dist;

		i = best-1;
	}

}

void insertion(tsp_instance* tsp_in, int* visited_nodes)
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
	double* costs = (double*)calloc((size_t)tsp_in->num_nodes, sizeof(double));

	int count = 2;

	visited_nodes[0] = indices[0];
	visited_nodes[1] = indices[1];

	node1[0] = indices[0];
	node2[0] = indices[1];
	node1[1] = indices[1];
	node2[1] = indices[0];

	costs[0] = costs[1] = max_dist;

	double best_cost = max_dist * 2;

	if (tsp_in->integerDist)
		tsp_in->bestCostI = (int)(best_cost + CAST_PRECISION);
	else
		tsp_in->bestCostD = best_cost;

	int i_best;
	
	for (; count < tsp_in->num_nodes ; count++)
	{
		double best_cost_h = DBL_MAX;
		int k_best;

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
	
	for (i = 0; i < tsp_in->num_nodes; i++)
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

void min_extra_mileage(tsp_instance* tsp_in, int count, int* visited_nodes, int* node1, int* node2, double* costs, int* i_best, int* k_best, double* best_cost_h, double* best_cost )
{
#ifdef GRASP

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

//visited_nodes controllare che non venga sovrascritto in maniera sbagliata
void greedy_refinement(tsp_instance* tsp_in, int* visited_nodes, double* cost) //aggiungere costo soluzione attuale (ottimo locale)
//e se il costo computato dal refinement per la
//nuova soluzione è diverso rispetto a cost e maggiore rispetto a tsp_in->cost, allore param cost = questo costo
//altrimenti (se diverso rispetto a cost e minore rispetto a tsp_in->cost) allora param cost = tsp_in-> cost = questo costo
{
	int* succ = calloc(tsp_in->num_nodes, sizeof(int));

	//int n_comps = 1;

	//define_tour(tsp_in, tsp_in->sol, succ, comp, &n_comps);//sistemare, da visited_nodes a succ
	/*int h = 0;
	for (; h < tsp_in->num_nodes; h++)
		printf("visited_nodes[%d] = %d\n", h, visited_nodes[h]);
	printf("%s", LINE);
	*/
	succ_construction(visited_nodes, succ, tsp_in->num_nodes);

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

		int count_edge = 0;
		int j = 0;
		for (; j < tsp_in->num_nodes && count_edge < 2; j++)
		{
			if (j!=i && j!=succ[i] && succ[j]!=i && succ[j]!=succ[i])
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
					cost_j_h = (double) x1;
					cost_i_j = (double) x2;
					cost_k_h = (double) x3;
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
					count_edge++;
					/*
					tsp_in->sol[xpos(tsp_in, i, j)] = 1.0;
					tsp_in->sol[xpos(tsp_in, succ[i], succ[j])] = 1.0;
					tsp_in->sol[xpos(tsp_in, i, succ[i])] = 0.0;
					tsp_in->sol[xpos(tsp_in, j, succ[j])] = 0.0;
					*/
					printf("[%d, %d] [%d,%d] -->> [%d,%d] [%d,%d]\n ", i, succ[i], j, succ[j], i, j, succ[i], succ[j]);

					if (tsp_in->integerDist)
						tsp_in->bestCostI += delta;
					else
						tsp_in->bestCostD += delta;


					int k = succ[i];
					int next = succ[i];

					int count = 0;
					int* orientation = (int*) calloc((size_t) tsp_in->num_nodes, sizeof(int));

					while (next != j)
					{
						orientation[count++] = next;
						next = succ[next];
					}

					next = succ[k];
					int next_count=0;

					while (next != j)
					{
						int temp = succ[next];
						succ[next] = orientation[next_count++];
						next = temp;
					}

					succ[i] = j;
					succ[k] = succ[j];
					succ[j] = orientation[count - 1];
					
					/*
					int begin = 0;
					int count2 = 0;
					int r = 0;

					int node = -1;
					do
					{
						node = succ[r];
						r = node;
						visited_nodes[++count2] = node;  //inserirlo nell'algoritmo
					} while (node != begin);

						int* succ1 = calloc(tsp_in->num_nodes, sizeof(int));
						succ_construction(visited_nodes, succ1, tsp_in->num_nodes);

						int* comp1 = calloc(tsp_in->num_nodes, sizeof(int));
						int s;
						for (s = 0; s < tsp_in->num_nodes; s++)
							comp1[s] = 1;

						

						int n_comps = 1;
						if (tsp_in->plot)
						{
							plot(tsp_in, succ1, comp1, &n_comps);

							free(succ1);
							free(comp1);
						}
					*/

						//break;
				}
			}
		}
	}

	/*
	i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		printf("[%d,%d]\n", i, succ[i]);
	}
	*/

	visited_nodes[0] = 0;
	int begin = 0;
	int count = 0;
	i = 0;

	int node = -1;
	do
	{
		node = succ[i];
		i = node;
		visited_nodes[++count] = node;  //inserirlo nell'algoritmo
		printf("count: %d\n", visited_nodes[count - 1]);  
	} 
	while (node != begin);

	//se la soluzione è migliorata, rispetto a quella in input, modificala
}

void vns(tsp_instance* tsp_in, int* visited_nodes)
{
	//se la soluzione migliora, questa va inserita dentro visited nodes
	srand(time(NULL));
	
	int* kopt_visited_nodes = (int*) calloc((size_t)tsp_in->num_nodes, sizeof(int));

	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		kopt_visited_nodes[i] = visited_nodes[i];
	}

	int num_iterations=0;

	for(;num_iterations<10; num_iterations++)
	{
		int t = rand() % tsp_in->num_nodes;
		int k = 1;
		int max_k = tsp_in->num_nodes / 2;
	
		for (; k < max_k; k++)
		{
			int next_index = (t+k) % tsp_in->num_nodes;
			int temp = kopt_visited_nodes[next_index];
			kopt_visited_nodes[next_index] = kopt_visited_nodes[t];
			kopt_visited_nodes[t] = temp;

			double cost=0.0;
			i = 0;
			for (; i < tsp_in->num_nodes; i++)
			{
				if (tsp_in->integerDist)
				{
					int i1 = i % tsp_in->num_nodes;
					int i2 = (i + 1) % tsp_in->num_nodes;
					int x;
					dist(kopt_visited_nodes[i1],kopt_visited_nodes[i2], tsp_in, &x);
					cost += (double)x;
				}
				else
				{
					double x;
					int i1 = i % tsp_in->num_nodes;
					int i2 = (i + 1) % tsp_in->num_nodes;
					dist(kopt_visited_nodes[i1], kopt_visited_nodes[i2], tsp_in, &x);
					cost += x;
				}
			}
			
			//refinement

			//se peggiora la soluzione rispetto a quella attuale torno a calcolare t
			//altrimenti aggiorno visited_nodes con il valore di kopt_visited_nodes e aumento il conto delle soluzioni locali trovate

		}
	}

	//scandire visited nodes e aggiornare tsp_in->sol
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

void succ_construction(int* visited_nodes, int* succ, int num_nodes)//succ deve essere già allocato
{
	int i;
	for (i = 0; i < num_nodes; i++)
		succ[visited_nodes[i]] = visited_nodes[(i+1) % num_nodes];

}