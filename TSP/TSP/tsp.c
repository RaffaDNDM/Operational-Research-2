#include "tsp.h"
#include "input.h"
#include "cplex_solver.h"
#include <cplex.h>

int main(int argc, char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);

	manage_input(&tsp_in);
	
	return 0;
}

void solution(tsp_instance* tsp_in)
{
	tsp_in->sol = (int*)calloc(((size_t)tsp_in->num_nodes) +1, sizeof(int));

	cplex_solver(tsp_in);
	dealloc_inst(tsp_in);
}

void set_params_and_solve(tsp_instance* tsp_in)
{
	parse_file(tsp_in);

	#ifndef METAHEURISTIC
		solution(tsp_in);
	#else
		int node_lim[] = { 0, 1, 2, 3 };
		int sol_lim[] = { 1, 2, 3 };
		double gap[] = { 0.1, 0.01, 0.001, 0.0001 };
		int seed[] = { 500, 1500, 2000, 2500 };
		
		#ifdef ALL_PARAM_COMBINATIONS
		
			double best_time = CPX_INFBOUND;
			int h = 0;
			for (; h < (sizeof(seed) / sizeof(int)); h++)
			{
				if (tsp_in->alg == 1)
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
								solution(tsp_in);

								if (tsp_in->execution_time < best_time)
									best_time = tsp_in->execution_time;
							}
						}
					}
				}
				else
				{
					tsp_in->seed = seed[h];
					solution(tsp_in);

					if (tsp_in->execution_time < best_time)
						best_time = tsp_in->execution_time;
				}
			}

			tsp_in->execution_time = best_time;
		#else
			tsp_in->node_lim = node_lim[0];
			tsp_in->sol_lim = sol_lim[0];
			tsp_in->eps_gap = gap[0];
			tsp_in->seed = seed[0];
			solution(tsp_in);
		#endif
	#endif
}

void manage_input(tsp_instance* tsp_in)
{
	char* name_algs[] = { ALG1, ALG2, ALG3, ALG4 };
	FILE* perf_data = NULL;

	if (strncmp(tsp_in->dir, "NULL", 4) != 0)
	{
		FILE* pipe = _popen(DIR_LIST_PY, "w");
		fprintf(pipe, "%s\n%d\n", tsp_in->dir, tsp_in->size);
		_pclose(pipe);
		
		//printf("Read File \n");
		//char c=getchar();
		FILE* instances = fopen("instances.txt", "r");

		#ifdef PERF_PROF_ON
			char title[40];
			int num_algs = 0;
			perf_data = fopen("perf_data.csv", "w");

			int i = 0;
			for (; i < NUM_ALGS; i++)
			{
				if (tsp_in->which_alg[i])
				{
					num_algs++;
				}
			}

			fprintf(perf_data, "%d", num_algs);
			i = 0;
			for (; i < NUM_ALGS; i++)
			{
				if (tsp_in->which_alg[i])
				{
					fprintf(perf_data, ", %s", name_algs[i]);
				}
			}
			fprintf(perf_data, "\n");
		#endif

		char in_file[LINE_SIZE];
		while (fgets(in_file, LINE_SIZE, instances) != NULL)
		{
			in_file[strlen(in_file) - 1] = 0;
			strcpy(tsp_in->input, in_file);
			fprintf(perf_data, "%s", tsp_in->input);

			if (tsp_in->alg == 0)
			{
				int n = 0;
				for (; n < NUM_ALGS; n++)
				{
					if (tsp_in->which_alg[n])
					{
						tsp_in->alg = n + 1;
						set_params_and_solve(tsp_in);

						#ifdef PERF_PROF_ON
						/*
						if (tsp_in.integerDist)
							fprintf(perf_data, ", %d", tsp_in->bestCostI);
						else
							fprintf(perf_data, ", %lf", tsp_in->bestCostI);
						*/
					fprintf(perf_data, ", %lf", tsp_in->execution_time);
					#endif
				}
			}

			#ifdef PERF_PROF_ON
				fprintf(perf_data, "\n");
			#endif

			tsp_in->alg = 0;
			
			}
			else
				set_params_and_solve(tsp_in);
		}

		#ifdef PERF_PROF_ON
			fclose(perf_data);
		#endif
	}
	else
	{
		if (tsp_in->alg == 0)
		{
			int n = 0;
			for (; n < NUM_ALGS; n++)
			{
				if (tsp_in->which_alg[n])
				{
					tsp_in->alg = n + 1;
					set_params_and_solve(tsp_in);
				}
			}
		}
		else
			set_params_and_solve(tsp_in);
	}
	/*
	FILE* pipe = _popen(PERF_PROF_PY, "w");
	_pclose(pipe);
	*/
}