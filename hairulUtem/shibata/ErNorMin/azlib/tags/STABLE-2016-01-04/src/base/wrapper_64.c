/*********************************************************************
 * wrapper_64.c
 *
 * Copyright (C) 2005 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: wrapper_64.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include "wrapper_64.h"
#include "rc.h"
#ifdef WIN32
#include <io.h>
#endif /* WIN32 */


WRAP64_INT
wrap64_ftell (FILE *fp)
{
	WRAP64_INT ret;

	RC_TRY( wrap64_fseek(fp, (WRAP64_INT)0, SEEK_CUR) );
	/* fflush(fp);*/

#ifdef WIN32
	ret = (WRAP64_INT)_telli64(_fileno(fp));
#else  /* WIN32 */
  #ifdef  _LARGEFILE_SOURCE
	ret = (WRAP64_INT)ftello(fp);
  #else
	ret = (WRAP64_INT)ftell(fp);
  #endif
#endif  /* WIN32 */

	return(ret);
}


RC
wrap64_fseek (FILE *fp, WRAP64_INT offset, int origin)
{
	if(fflush(fp) != 0) return(FLUSH_ERROR_RC);
#ifdef WIN32
	{   
		WRAP64_INT ret;

		if(origin == SEEK_SET) rewind(fp); /* なぜかrewind()が必要（バグ？） */
		ret = _lseeki64(_fileno(fp), (__int64)offset, origin);
		if(ret < (WRAP64_INT)(0L)) return(SEEK_ERROR_RC);
	}
#else  /* WIN32 */
  #ifdef  _LARGEFILE_SOURCE
	if(fseeko(fp, (off_t)offset, origin) != 0) return(SEEK_ERROR_RC);
  #else
	if(fseek(fp, (long)offset, origin) != 0) return(SEEK_ERROR_RC);
  #endif
#endif  /* WIN32 */

	return(NORMAL_RC);
}


