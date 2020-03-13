#include "tsp.h"
#include "input.h"

int main(char** argv, int argc)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);
	parse_file(&tsp_in);
	
	return 0;
}