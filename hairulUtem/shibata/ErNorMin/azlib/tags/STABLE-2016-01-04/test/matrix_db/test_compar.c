#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"
#include "matrix_db.h"

int main (int argc, char **argv)
{
	fprintf(stderr, "row1, col1, row2, col2 \n");
	fprintf(stderr, "1,2,2,3 = %d\n", compar_row_col(1,2,2,3));
	fprintf(stderr, "1,3,2,2 = %d\n", compar_row_col(1,3,2,2));
	fprintf(stderr, "2,1,2,1 = %d\n", compar_row_col(2,1,2,1));
	fprintf(stderr, "2,1,3,1 = %d\n", compar_row_col(2,1,3,1));
	fprintf(stderr, "3,1,2,1 = %d\n", compar_row_col(3,1,2,1));
	fprintf(stderr, "2,2,1,3 = %d\n", compar_row_col(2,2,1,3));
	fprintf(stderr, "2,1,2,2 = %d\n", compar_row_col(2,1,2,2));
	fprintf(stderr, "2,2,2,1 = %d\n", compar_row_col(2,2,2,1));

	return(EXIT_SUCCESS);
}


