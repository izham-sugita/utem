#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "base.h"
#include "mathematics.h"

int main (int argc, char **argv)
{
	FILE *fp;

	if(argc < 2){
		fprintf(stderr, "Usage : %s [r:large_file]\n", argv[0]);
		fprintf(stderr, " Please make a large file(> 2GB)\n");
		fprintf(stderr, "  ex) $cat /dev/zero >> file (stop CTRL+c)\n");
		return(EXIT_FAILURE);
	}

	/* Log Print Level */
	RC_TRY_MAIN( set_log_level(5) );

	/* データモデル チェック */
	print_data_model(stderr);

	/* 2GB以上ファイルのファイルが開けるか確認 */
	RC_TRY_MAIN( log_printf(1, "Opne Large File...") );
	RC_TRY_MAIN( rc_fopen(argv[1], "r", &fp) );
	RC_TRY_MAIN( log_printf(1, "OK!!\n") );
	RC_TRY_MAIN( rc_fclose(fp) );
	
	return(EXIT_SUCCESS);
}

