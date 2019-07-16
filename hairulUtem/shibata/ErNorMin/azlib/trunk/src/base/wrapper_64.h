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

/* $Id: wrapper_64.h 842 2006-12-04 06:21:20Z sasaoka $ */

#ifndef WRAPPER_64_H
#define WRAPPER_64_H

#include <stdio.h>
#include "rc.h"

/* 64bit 整数型の定義
 * C99 対応コンパイラなら long long を使えば良いのだが，時機尚早．
 * Win32環境以外の 64bit OS(主にUNIX) では long を 64bit と仮定し，
 * 32bitのLinux等は _LARGEFILE_SOURCE および _FILE_OFFSET_BITS=64 が
 * 定義されていれば off_t を 64bit 整数として使用
 */
#ifdef WIN32
typedef __int64 WRAP64_INT;
#else /* WIN32 */
 #ifdef _LARGEFILE_SOURCE
   typedef off_t WRAP64_INT; /* 32bit Linux用 */
 #else
   typedef long WRAP64_INT; /* 64bit UNIX用 */
 #endif 
#endif  /* WIN32 */

/* 32bit/64bit 切り替え可能な整数の定義 */
#ifdef INTEGER_SIZE64
#ifdef WIN32
typedef __int64  WRAP32_64_INT;
#else  /* WIN32 */
typedef long     WRAP32_64_INT;
#endif /* WIN32 */
#else  /* INTEGER_SIZE64 */
#ifdef WIN32
typedef int      WRAP32_64_INT;
#else  /* WIN32 */
typedef int      WRAP32_64_INT;
#endif /* WIN32 */
#endif /* INTEGER_SIZE64 */

/* wrapper_64.c */
WRAP64_INT wrap64_ftell(FILE *fp);
RC wrap64_fseek(FILE *fp, WRAP64_INT offset, int origin);


#endif  /* WRAPPER_64 */

