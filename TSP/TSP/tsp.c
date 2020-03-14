#include "tsp.h"
#include "input.h"
#include <assert.h>
#include <stdio.h>

int main(int argc,char** argv)
{
	tsp_instance tsp_in;
	parse_cmd(argv, argc, &tsp_in);
	parse_file(&tsp_in);

	dealloc_inst(&tsp_in);
	
	return 0;
}