/*********************************************************************
 * mtrix_db_plus.c
 *
 * Copyright (C) 2006 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Thanks to following contributors.
 *   <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id$ */

#include <stdio.h>
#include <math.h>
#ifndef WIN32
#include <unistd.h>
#endif  /* WIN32 */
#include <string.h>
#include <limits.h>
#include "wrapper_64.h"
#include "scratch_io.h"
#include "matrix_db_plus.h"
#include "mathematics.h"
#include "rc.h"
#include "log_printf.h"
#include "memory_manager.h"
#include "macros.h"


void
print_matrix_db_size (FILE *fp)
{
	fprintf(fp, " MATRIX_DB_INDEX_SIZE   = %d\n", MATRIX_DB_INDEX_SIZE);
	fprintf(fp, " MATRIX_DB_INDEX_HEADER = %d\n", MATRIX_DB_INDEX_HEADER);
#ifdef WIN32
	fprintf(fp, " MATRIX_DB_INDEX        = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_INDEX) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, " MATRIX_DB_INDEX        = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_INDEX) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, " MATRIX_DB_INDEX        = %ld [Byte]\n",
	                                       sizeof(MATRIX_DB_INDEX) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */
	fprintf(fp, " MATRIX_DB_DATA_SIZE    = %d\n", MATRIX_DB_DATA_SIZE);
	fprintf(fp, " MATRIX_DB_DATA_HEADER  = %d\n", MATRIX_DB_DATA_HEADER);
#ifdef WIN32
	fprintf(fp, " MATRIX_DB_DATA         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_DATA) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, " MATRIX_DB_DATA         = %d [Byte]\n",
	                                       sizeof(MATRIX_DB_DATA) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, " MATRIX_DB_DATA         = %ld [Byte]\n",
	                                       sizeof(MATRIX_DB_DATA) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */
}

