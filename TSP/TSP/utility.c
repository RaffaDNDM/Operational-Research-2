#include "utility.h"
#include <math.h>

void dist(int node1, int node2, tsp_instance* tsp_in, void* dist)
{
	double x_dist = tsp_in->x_coords[node1] - tsp_in->x_coords[node2];
	double y_dist = tsp_in->y_coords[node1] - tsp_in->y_coords[node2];

	if (tsp_in->integerDist)
	{
		int x_distI = (int)(x_dist + CAST_PRECISION);
		int y_distI = (int)(y_dist + CAST_PRECISION);

		int* distI = (int*)dist;
		*distI = (int)(sqrt((double)x_distI * (double)x_distI + (double)y_distI * (double)y_distI) + CAST_PRECISION);

		return;
	}

	double* distD = (double*)dist;
	*distD = sqrt(x_dist * x_dist + y_dist * y_dist);
}

void print_cost(tsp_instance* tsp_in)
{
	printf(LINE);
	printf("The solution cost is:\n");
	if (tsp_in->integerDist)
		printf("int cost %d\n", tsp_in->bestCostI);
	else
		printf("double cost %10.31f\n", tsp_in->bestCostD);
	printf(LINE);
}