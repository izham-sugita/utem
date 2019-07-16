/*********************************************************************
 * data_model.c
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: data_model.c 858 2006-12-16 03:32:08Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rc.h"
#include "wrapper_64.h"


void
print_data_model (FILE *fp)
{

	fprintf(fp, "              bit byte\n");
#ifdef WIN32
	fprintf(fp, "char        = %3d  %3d\n"
	            "short       = %3d  %3d\n"
	            "int         = %3d  %3d\n"
	            "long        = %3d  %3d\n"
	            "float       = %3d  %3d\n"
	            "double      = %3d  %3d\n"
	            "long double = %3d  %3d\n"
	            "void *      = %3d  %3d\n"
	            "size_t      = %3d  %3d\n"
	            "time_t      = %3d  %3d\n"
	            "WRAP64_INT  = %3d  %3d\n", 
	            8*sizeof(char),        sizeof(char),
	            8*sizeof(short),       sizeof(short),
	            8*sizeof(int),         sizeof(int),
	            8*sizeof(long),        sizeof(long),
	            8*sizeof(float),       sizeof(float),
	            8*sizeof(double),      sizeof(double),
	            8*sizeof(long double), sizeof(long double),
	            8*sizeof(void *),      sizeof(void *),
	            8*sizeof(size_t),      sizeof(size_t),
	            8*sizeof(time_t),      sizeof(time_t),
	            8*sizeof(WRAP64_INT),  sizeof(WRAP64_INT) );
#else /* WIN32 */
  #ifdef _LARGEFILE_SOURCE
	fprintf(fp, "char        = %3d  %3d\n"
	            "short       = %3d  %3d\n"
	            "int         = %3d  %3d\n"
	            "long        = %3d  %3d\n"
    #ifdef USE_ISOC99 
	            "long long   = %3d  %3d\n"
    #endif /* USE_ISOC99 */
	            "float       = %3d  %3d\n"
	            "double      = %3d  %3d\n"
	            "long double = %3d  %3d\n"
	            "void *      = %3d  %3d\n"
	            "size_t      = %3d  %3d\n"
	            "time_t      = %3d  %3d\n"
	            "WRAP64_INT  = %3d  %3d\n",
	            8*sizeof(char),        sizeof(char),
	            8*sizeof(short),       sizeof(short),
	            8*sizeof(int),         sizeof(int),
	            8*sizeof(long),        sizeof(long),
    #ifdef USE_ISOC99 
	            8*sizeof(long long),   sizeof(long long),
    #endif /* USE_ISOC99 */
	            8*sizeof(float),       sizeof(float),
	            8*sizeof(double),      sizeof(double),
	            8*sizeof(long double), sizeof(long double),
	            8*sizeof(void *),      sizeof(void *),
	            8*sizeof(size_t),      sizeof(size_t),
	            8*sizeof(time_t),      sizeof(time_t),
	            8*sizeof(WRAP64_INT),  sizeof(WRAP64_INT) );
  #else /* _LARGEFILE_SOURCE */
	fprintf(fp, "char        = %3ld  %3ld\n"
	            "short       = %3ld  %3ld\n"
	            "int         = %3ld  %3ld\n"
	            "long        = %3ld  %3ld\n"
    #ifdef USE_ISOC99 
	            "long long   = %3ld  %3ld\n"
    #endif /* USE_ISOC99 */
	            "float       = %3ld  %3ld\n"
	            "double      = %3ld  %3ld\n"
	            "long double = %3ld  %3ld\n"
	            "void *      = %3ld  %3ld\n"
	            "size_t      = %3ld  %3ld\n"
	            "time_t      = %3ld  %3ld\n"
	            "WRAP64_INT  = %3ld  %3ld\n",
	            8*sizeof(char),        sizeof(char),
	            8*sizeof(short),       sizeof(short),
	            8*sizeof(int),         sizeof(int),
	            8*sizeof(long),        sizeof(long),
    #ifdef USE_ISOC99 
	            8*sizeof(long long),   sizeof(long long),
    #endif /* USE_ISOC99 */
	            8*sizeof(float),       sizeof(float),
	            8*sizeof(double),      sizeof(double),
	            8*sizeof(long double), sizeof(long double),
	            8*sizeof(void *),      sizeof(void *),
	            8*sizeof(size_t),      sizeof(size_t),
	            8*sizeof(time_t),      sizeof(time_t),
	            8*sizeof(WRAP64_INT),  sizeof(WRAP64_INT) );
  #endif /* _LARGEFILE_SOURCE */
#endif  /* WIN32 */

}

