#include "tsp.h"
#include "input.h"
#include "algorithms.h"
#include <cplex.h>

int main(int argc,char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);
	parse_file(&tsp_in);
	solution(&tsp_in);		

	dealloc_inst(&tsp_in);
	
	return 0;
}

void solution(tsp_instance* tsp_in)
{
	tsp_in->sol = (int*)calloc(((size_t)tsp_in->num_nodes) +1, sizeof(int));

	//default behaviour 
	if (tsp_in->alg == 1)
		default_alg(tsp_in);

	if (tsp_in->alg == 2)
		cplex_solver(tsp_in);
}


void dist(int node1, int node2, tsp_instance* tsp_in, void* dist)
{
	double x_dist = tsp_in->x_coords[node1] - tsp_in->x_coords[node2];
	double y_dist = tsp_in->y_coords[node1] - tsp_in->y_coords[node2];
	
	if (tsp_in->integerDist)
	{	
		int x_distI = (int) (x_dist + CAST_PRECISION);
		int y_distI = (int) (y_dist + CAST_PRECISION);
		
		int* distI = (int*) dist;
		*distI =(int)( sqrt((double)x_distI *(double) x_distI +(double)y_distI *(double) y_distI) + CAST_PRECISION);

		return;
	}

	double* distD = (double*) dist;
	*distD = sqrt(x_dist*x_dist + y_dist*y_dist);
}