/*********************************************************************
 * wrapper_64.h
 *
 * Copyright (C) 2006 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: wrapper_64.h 422 2005-08-19 09:52:15Z sasaoka $ */

#ifndef WRAPPER_64_H
#define WRAPPER_64_H

#include <stdio.h>
#include "rc.h"

/* 64bit 整数型の定義 */
/* C99 対応コンパイラなら long long を使えば良いのだが... */
/* Win32環境以外は 64bit 環境と仮定 */
#ifdef WIN32
typedef __int64 WRAP64_INT;
#else /* WIN32 */
 #ifdef _LARGEFILE_SOURCE
   typedef off_t WRAP64_INT;
 #else
   typedef long WRAP64_INT;
 #endif 
#endif  /* WIN32 */

/* wrapper_64.c */
WRAP64_INT wrap64_ftell(FILE *fp);
RC wrap64_fseek(FILE *fp, WRAP64_INT offset, int origin);


#endif  /* WRAPPER_64 */

