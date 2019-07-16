#include <stdio.h>
#include <stdlib.h>
#include "base.h"
#include "mathematics.h"


int main(void);

int main(void)
{
	int ii1,ii2;
	char *c = NULL;

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );
	RC_TRY_MAIN( set_log_file(0, stderr) );

	print_data_model(stderr);

	ii1 = 1;
	RC_TRY_MAIN( log_printf(5, "MEMORY ALLOCATE TEST START\n"
	                           "STOP : [Ctrl + C]\n") );
	while(1){
		c = (char *)realloc(c, (size_t)ii1*MEGA_BYTE*sizeof(char));
		if(c == (char *)NULL){
			RC_TRY_MAIN( error_printf(9001, ii1) );
			return(EXIT_FAILURE);
		}

		for(ii2=0;ii2<MEGA_BYTE;ii2++){
			c[((size_t)ii1-1)*MEGA_BYTE+ii2] = 'a';
		}
		RC_TRY_MAIN( log_printf(1, "Now %4d MByte allocated!!\r", ii1) );
		ii1++;
	}

	return(EXIT_SUCCESS);
}

