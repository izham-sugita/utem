#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "mathematics.h"
#include "fem.h"
#include "base.h"



int main(int argc, char **argv);

int main(int argc, char **argv)
{
	FILE *fp;
	NODE_ARRAY node;
	ELEMENT_ARRAY element;
	RIGID_ELEMENT_ARRAY rigid;
	MATERIAL_PROP_ARRAY material;
	PHYSICAL_PROP_ARRAY physical;
	DISP_ARRAY mode[10];
	BC_ARRAY rest;
	double evals[10];
	int e_num;
	int ii1;


	if(argc < 2){
		fprintf(stderr, "Usage : %s [r:model(.bdf)]\n", argv[0]);
		return(EXIT_FAILURE);
	}

	/* MM, Log */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( mm_init((size_t)600*MEGA_BYTE) );
	RC_TRY_MAIN( set_log_file(0, stderr) );
	

	/* Model Input */
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node) );
	RC_TRY_MAIN( nst_input_element(fp, &element) );
	RC_TRY_MAIN( nst_input_material_prop(fp, &material) );
	RC_TRY_MAIN( nst_input_physical_prop(fp, &physical) );
	RC_TRY_MAIN( nst_input_rigid_element(fp, &rigid) );
	RC_TRY_MAIN( nst_input_restraint(fp, &rest) );
	RC_TRY_MAIN( rc_fclose(fp) );

//	RC_TRY_MAIN( cut_element(&node, &element) );
//	RC_TRY_MAIN( dummy_renumber(&node) );
	RC_TRY_MAIN( fem_renumber0(&node, element, rigid, 0) );

	RC_TRY_MAIN( log_printf(5, "node = %d\n", node.size) );
	RC_TRY_MAIN(fem_eigen_block_solve_sky3(node, element, material,
	                                       physical, rest, 1000.0, 
                                           7, 20, 10, &e_num, evals,
	                                       mode) );
	for(ii1=0; ii1<e_num; ii1++){
		RC_TRY_MAIN( log_printf(5, "e[%d] = %e\n", ii1,evals[ii1]) );
	}

	

	RC_TRY_MAIN( free_node_array(&node) );
	RC_TRY_MAIN( free_element_array(&element) );
	RC_TRY_MAIN( free_material_prop_array(&material) );
	RC_TRY_MAIN( free_physical_prop_array(&physical) );

	RC_TRY_MAIN( mm_terminate() );

	return(EXIT_SUCCESS);
}
