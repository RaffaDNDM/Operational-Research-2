#include "tsp.h"
#include "input.h"
#include "cplex_solver.h"
#include "default_alg.h"
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

	switch (tsp_in->alg)
	{
		//default behaviour 
		case 1:
			default_alg(tsp_in);
			break;

		case 2:
			cplex_solver(tsp_in);
			break;
	}
	
}