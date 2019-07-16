#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "base.h"
#include "mathematics.h"
#include "fem.h"
#include "matrix_db.h"

int
main(int argc, char **argv)
{
	int ii1;
	FILE *fp;
	ELEMENT_ARRAY elem1, elem2, surf1, surf2;
	NODE_ARRAY node1, node2;
	int size1, size2;
	int *index_array1, *index_array2;

	if(argc < 3){
		fprintf(stderr, "Usage : %s [r:model1(.NAS)] [r:model2(.NAS)]\n",
		        argv[0]);
		return(EXIT_FAILURE);
	}

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( set_log_file(0, stderr) )

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)1024*MEGA_BYTE) );

	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node1) );
	RC_TRY_MAIN( nst_input_element(fp, &elem1) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( dummy_renumber(&node1) );
	log_printf(5, "node1.size = %d\n", node1.size);
	log_printf(5, "elem1.size = %d\n", elem1.size);
	
	RC_TRY_MAIN( rc_fopen(argv[2], "r", &fp) );
	RC_TRY_MAIN( nst_input_node(fp, &node2) );
	RC_TRY_MAIN( nst_input_element(fp, &elem2) );
	RC_TRY_MAIN( rc_fclose(fp) );
	RC_TRY_MAIN( dummy_renumber(&node2) );
	log_printf(5, "node2.size = %d\n", node2.size);
	log_printf(5, "elem2.size = %d\n", elem2.size);
#if 1
	{
		for(ii1=0; ii1<node2.size; ii1++){
			node2.array[ii1].p.y -= 0.001f; /* キューブ */
			//node2.array[ii1].p.x += 200.0f; /* Knuckle */
		}
	}
#endif

	DEBUG_PRINT;
	RC_TRY_MAIN( extract_surface(elem1, &surf1) );
	RC_TRY_MAIN( extract_surface(elem2, &surf2) );
	log_printf(5, "surf1.size = %d\n", surf1.size);
	log_printf(5, "surf2.size = %d\n", surf2.size);

	DEBUG_PRINT;
	RC_NULL_CHK( index_array1 = mm_alloc(surf1.size*sizeof(int)) );
	RC_NULL_CHK( index_array2 = mm_alloc(surf2.size*sizeof(int)) );

	DEBUG_PRINT;
	timer_start(0);
	RC_TRY_MAIN( extract_geom_contact_surface(surf1, node1, surf2, node2,
	                                          &size1, index_array1,
	                                          &size2, index_array2) );
	RC_TRY_MAIN( log_printf(5, "time : %d\n", timer_end(0)) );
	RC_TRY_MAIN( log_printf(5, "size1 =  %d\n", size1) );
	RC_TRY_MAIN( log_printf(5, "size2 =  %d\n", size2) );
	DEBUG_PRINT;

	RC_TRY_MAIN( mm_free((void *)index_array1) );
	RC_TRY_MAIN( mm_free((void *)index_array2) );
	RC_TRY_MAIN( free_node_array(&node1) );
	RC_TRY_MAIN( free_node_array(&node2) );
	RC_TRY_MAIN( free_element_array(&elem1) );
	RC_TRY_MAIN( free_element_array(&elem2) );

	/* Memory Manager */
	RC_TRY_MAIN( mm_terminate() );

	return(EXIT_SUCCESS);
}

