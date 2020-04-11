#include "tsp.h"
#include "input.h"
#include "cplex_solver.h"
#include "default_alg.h"
#include <cplex.h>

int main(int argc,char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);
	/* if (tsp.dir=! null o strings vuota ) allora c'è una cartella
	{
		leggi(tsp.dir) dovrebbe restituire array di stringhe , nomi dei file delle istanze 
		while
		{
			tsp_in.input=istanza letta array_istanze[i]
			parse_file(&tsp_in);
			solution(&tsp_in);
		}
	}
	else 
	{
		parse_file(&tsp_in); c'era solo un file
		solution(&tsp_in);
	}
	*/





	parse_file(&tsp_in);
	
	//lettura file da input (funzioni per qualsiasi sistema operativo), ripetere solution su ogni istanza
	//array da inserire in cicli interni che provi tutte le combinazioni
	//metterli dentro tsp_in
	int node_lim[] = { 0, 1, 2, 3 };
	int sol_lim[] = { 1, 2, 3 };
	double gap[] = { 0.1, 0.01, 0.001, 0.0001 };
	//aggiungere array seed
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