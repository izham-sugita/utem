#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "rc.h"
#include "math_utl.h"
#include "macros.h"


#define SIZE (10000)

int main(void);

int main(void)
{
	int ii1;

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( set_log_file(0, stderr) )

	fprintf(stderr, "Calculating ");
	for(ii1=0; ii1<SIZE; ii1++) {
		RC_TRY_MAIN( progress_bar_printf(5, (double)(ii1+1)/SIZE,
		                                  " %3d/%3d", (ii1+1), SIZE));

	}
	fprintf(stderr, "\n");

	return(EXIT_SUCCESS);
}

