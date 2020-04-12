#include "tsp.h"
#include "input.h"
#include "cplex_solver.h"
#include "default_alg.h"
#include <cplex.h>
//#include <windows.h>

int main(int argc, char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);

	if (tsp_in.dir != NULL)
	{	
		FILE* pipe = _popen(DIR_LIST_PY, "w");
		fprintf(pipe, "%s", tsp_in.dir);
		_pclose(pipe);

		FILE* instances = fopen("instances.txt", "r");
		char* in_file = malloc(sizeof(char)*LINE_SIZE);
		while (fgets(in_file, LINE_SIZE, instances) != NULL)
		{
			in_file[strlen(in_file) - 1] = 0;
			strcpy(tsp_in.input, in_file);
 
			if (tsp_in.model == 0)
			{
				int n = 0;
				for (; n < NUM_MODELS; n++)
				{
					tsp_in.model = n + 1;
					set_params_and_solve(&tsp_in);
				}
			}
			else
				set_params_and_solve(&tsp_in);

			dealloc_inst(&tsp_in);
		
		}

	}
	else
	{
		if (tsp_in.model == 0)
		{
			int n = 0;
			for (; n < NUM_MODELS; n++)
			{
				tsp_in.model = n + 1;
				set_params_and_solve(&tsp_in);
			}
		}
		else
			set_params_and_solve(&tsp_in);
	}
	

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

void set_params_and_solve(tsp_instance* tsp_in)
{
	int node_lim[] = { 0, 1, 2, 3 };
	int sol_lim[] = { 1, 2, 3 };
	double gap[] = { 0.1, 0.01, 0.001, 0.0001 };
	int seed[] = { 500, 1500, 2000, 2500 };

	int h = 0;
	for (; h < (sizeof(seed) / sizeof(int)); h++)
	{
		if (tsp_in->alg == 2 && tsp_in->model == 1)
		{
			int i = 0;
			for (; i < (sizeof(node_lim) / sizeof(int)); i++)
			{
				int j = 0;
				for (; j < (sizeof(sol_lim) / sizeof(int)); j++)
				{
					int k = 0;
					for (; k < (sizeof(gap) / sizeof(double)); k++)
					{
						tsp_in->node_lim = node_lim[i];
						tsp_in->sol_lim = sol_lim[j];
						tsp_in->eps_gap = gap[k];
						tsp_in->seed = seed[h];
						parse_file(tsp_in);
						solution(tsp_in);
					}
				}
			}
		}
		else
		{
			tsp_in->seed = seed[h];
			parse_file(tsp_in);
			solution(tsp_in);
		}

	}
}