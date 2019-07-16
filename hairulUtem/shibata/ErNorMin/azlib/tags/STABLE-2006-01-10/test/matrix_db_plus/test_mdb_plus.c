#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "base.h"
#include "mathematics.h"
#include "matrix_db_plus.h"

int
main (int argc, char **argv)
{
	MATRIX_DB_INFO mdbi;
	MATRIX_DB_INDEX index;
	MATRIX_DB_DATA data;
	MATRIX_DB_SCRATCH scratch;
	MATRIX_DB_CACHE cache;

	/*
	if(argc < 2){
		fprintf(stderr, "Usage : %s [scratch directory(ex:/tmp)]\n", argv[0]);
		return(EXIT_FAILURE);
	}
	*/

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );

	/* Memory Manager */
	RC_TRY_MAIN( mm_init((size_t)256*MEGA_BYTE) );

	print_matrix_db_size(stderr);

	/* Memory Manager */
	RC_TRY( mm_terminate() );

	return(EXIT_SUCCESS);
}

