#include "tsp.h"
#include "input.h"
#include "algorithms.h"

int main(int argc,char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);
	parse_file(&tsp_in);
	solution(&tsp_in);
	
	if(tsp_in.plot)
		plot_solution(&tsp_in);

	dealloc_inst(&tsp_in);
	
	return 0;
}

void solution(tsp_instance* tsp_in)
{
	tsp_in->sol = (int*)calloc(tsp_in->num_nodes, sizeof(int));

	//default behaviour 
	if (tsp_in->alg == 1)
		default_alg(tsp_in);

	evaluate_sol(tsp_in);

	
	printf(LINE);
	printf("The solution cost is:\n");
	
	if (tsp_in->integerDist)
		printf("int cost %d\n", tsp_in->bestCostI);
	else
		printf("double cost %10.31f\n",tsp_in->bestCostD);
	
	printf(LINE);
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
		for (; i < tsp_in->num_nodes; i++)
		{
			int distI;

			dist(i - 1, i, tsp_in, &distI, 1);
			costI += distI;
		}
		tsp_in->bestCostI = costI;
	}
	else
	{
		for (; i < tsp_in->num_nodes; i++)
		{
			double distD;

			dist(i - 1, i, tsp_in, &distD, 0);
			costD += distD;
		}
		tsp_in->bestCostD = costD;
	}
}

void dist(int node1, int node2, tsp_instance* tsp_in, void* dist,  int costInt)
{
	double x_dist = tsp_in->x_coords[node1] - tsp_in->x_coords[node2];
	double y_dist = tsp_in->y_coords[node1] - tsp_in->y_coords[node2];
	
	if (costInt)
	{	
		int x_distI = (int) (x_dist + CAST_PRECISION);
		int y_distI = (int) (y_dist + CAST_PRECISION);
		
		int* distI = (int*) dist;
		*distI = sqrt(x_distI * x_distI + y_distI * y_distI) + CAST_PRECISION;

		return;
	}

	double* distD = (double*) dist;
	*distD = sqrt(x_dist*x_dist + y_dist*y_dist);
}