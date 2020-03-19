#include "algorithms.h"

void default_alg(tsp_instance* tsp_in)
{
	int i = 0;
	for (; i < tsp_in->num_nodes; i++)
	{
		tsp_in->sol[i] = i;
	}
}