/*********************************************************************
 * rc.c
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: rc.c 1043 2016-02-28 23:58:30Z hayashi $ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "rc.h"

#define FUNC_SIZE 10

static size_t cal_max_nmemb(size_t size);

/* 引数を 4つ持つ関数のポインタの配列 */
static void (*RcErrorPrintFunc[FUNC_SIZE])(const char *, int, const char *, RC)
  = {rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default,
     rc_error_print_default};

static int StackPtr = 0;   /* RcErrorPrintFunc[] 管理用のスタックポインタ */

/* 標準入出力関数 fopen() をRC型で定義した関数．
 * ファイルパス path で示されるファイルを，指定モード mode で開く．
 * 開いたファイルのポインタを fp に格納する．
 */
RC
rc_fopen (const char *path, const char *mode, FILE **fp)
{
	errno = 0;
	*fp = fopen(path, mode);

	/* オープンエラーの場合 */
	if(*fp == NULL){
		RC rc = OPEN_ERROR_RC;
		rc_error_print(__FILE__, __LINE__, path, rc);
		fprintf(stderr, "%s\n", strerror(errno));

		/* 大域変数 errno の表示 (5: verbose log mode) */
		/* 以下 LFS(Large File Summit)仕様 */
#ifdef _LAERGEFILE_SOURCE
		switch(errno)
		{
			case EFBIG:
				RC_TRY( error_printf(1004, "fopen(): errno = EFBIG") );
				break;
			case EOVERFLOW:
				RC_TRY( error_printf(1004, "fopen(): errno = EOVERFLOW") );
				break;
			default:
				RC_TRY( error_printf(1004, "fopen(): errno = unknown") );
		}
#endif
		return(rc);
	}

	return(NORMAL_RC);
}

/* 標準入出力関数 fclose() をRC型で定義した関数．
 * ファイルポインタ fp で示されるファイルを閉じる．
 */
RC
rc_fclose (FILE *fp)
{
	int ret;

	ret = fclose(fp);
	if(ret != 0) return(CLOSE_ERROR_RC);

	return(NORMAL_RC);
}

/* ファイルパス path で示されるファイルを作成する．
 */
RC
rc_touch (const char *path)
{
	FILE *fp;

	RC_TRY( rc_fopen(path, "a", &fp) );
	RC_TRY( rc_fclose(fp) );

	return(NORMAL_RC);
}


/* 出力先 stream を指定して，RCの列挙型変数名を出力する．
 */
void
print_rc(FILE *stream, RC rc)
{
	switch(rc){
	case NORMAL_RC:
		fprintf(stream,"NORMAL_RC\n");
		break;
	case END_RC:
		fprintf(stream,"END_RC\n");
		break;
	case IMPLEMENT_ERROR_RC:
		fprintf(stream,"IMPLEMENT_ERROR_RC\n");
		break;
	case READ_ERROR_RC:
		fprintf(stream,"READ_ERROR_RC\n");
		break;
	case WRITE_ERROR_RC:
		fprintf(stream,"WRITE_ERROR_RC\n");
		break;
	case CONVERT_ERROR_RC:
		fprintf(stream,"CONVERT_ERROR_RC\n");
		break;
	case NULL_ERROR_RC:
		fprintf(stream,"NULL_ERROR_RC\n");
		break;
	case ALLOC_ERROR_RC:
		fprintf(stream,"ALLOC_ERROR_RC\n");
		break;
	case OPEN_ERROR_RC:
		fprintf(stream,"OPEN_ERROR_RC\n");
		break;
	case CLOSE_ERROR_RC:
		fprintf(stream,"CLOSE_ERROR_RC\n");
		break;
	case SEEK_ERROR_RC:
		fprintf(stream,"SEEK_ERROR_RC\n");
		break;
	case FLUSH_ERROR_RC:
		fprintf(stream,"FLUSH_ERROR_RC\n");
		break;
	case ARG_ERROR_RC:
		fprintf(stream,"ARG_ERROR_RC\n");
		break;
	case CAL_ERROR_RC:
		fprintf(stream,"CAL_ERROR_RC\n");
		break;
	case OVERFLOW_ERROR_RC:
		fprintf(stream,"OVERFLOW_ERROR_RC\n");
		break;
	case UNDERFLOW_ERROR_RC:
		fprintf(stream,"UNDERFLOW_ERROR_RC\n");
		break;
	case ZERO_ERROR_RC:
		fprintf(stream,"ZERO_ERROR_RC\n");
		break;
	case NEGATIVE_ERROR_RC:
		fprintf(stream,"NEGATIVE_ERROR_RC\n");
		break;
	case SUBMIT_ERROR_RC:
		fprintf(stream,"SUBMIT_ERROR_RC\n");
		break;
	case OVERRUN_ERROR_RC:
		fprintf(stream,"OVERRUN_ERROR_RC\n");
		break;
	case IGNORE_ERROR_RC:
		fprintf(stream,"IGNORE_ERROR_RC\n");
		break;
	case UNKNOWN_ERROR_RC:
		fprintf(stream,"UNKNOWN_ERROR_RC\n");
		break;
	case SPECIAL_RC:
		fprintf(stream,"SPECIAL_RC\n");
		break;
	default:
		fprintf(stream,"???????\n");
	}
}


/* デフォルトのエラー出力関数．
 * ファイル名 file_name，行番号 line，関数名 rc_func，RC型 rc を受け取り，
 * それらを標準エラー出力に表示する．
 */
void
rc_error_print_default(const char *file_name, int line,
                       const char *rc_func, RC rc)
{
	fprintf(stderr, "[%s : %d] %s\n", file_name, line, rc_func);
	fprintf(stderr, "    ==> ");
	print_rc(stderr, rc);
}


/* デフォルトのエラー出力関数．
 * 何も表示しない．
 */
void
rc_error_print_silent(const char *file_name, int line,
                      const char *rc_func, RC rc)
{
	return;
}


/* エラー出力関数を func に変更する．
 */
void
rc_push_func(void (*func)(const char *, int, const char *, RC))
{
	StackPtr++;
	if(StackPtr >= FUNC_SIZE) StackPtr = FUNC_SIZE - 1;

	if(func == NULL) func = rc_error_print_default;
	RcErrorPrintFunc[StackPtr] = func;
}


/* エラー出力関数を元に戻す．
 */
void
rc_pop_func(void)
{
	StackPtr--;
	if(StackPtr < 0) StackPtr = 0;
}


/* エラー出力関数のドライバ．
 */
void
rc_error_print(const char *file_name, int line, const char *rc_func, RC rc)
{
	(RcErrorPrintFunc[StackPtr])(file_name, line, rc_func, rc);
}


static size_t
cal_max_nmemb(size_t size)
{
	if(size <= 1){
		return(1024*1024);
	}else if(size <= 2){
		return( 512*1024);
	}else if(size <= 4){
		return( 256*1024);
	}else if(size <= 8){
		return( 128*1024);
	}else if(size <= 128){
		return(  64*1024);
	}else if(size <= 1024){
		return(   8*1024);
	}

	return(1024);
}


/* fwrite() の RC 型ラッパー．
 */
RC
rc_fwrite(const void *ptr, size_t size, size_t nmemb, FILE *fp)
{
	RC_TRY( rc_fwrite_fn(ptr, size, nmemb, fp, NULL) );

	return(NORMAL_RC);
}


/* fwrite() の RC 型ラッパー．
 * エラー時にファイル名を表示する．
 */
RC
rc_fwrite_fn(const void *ptr, size_t size, size_t nmemb, FILE *fp,
             const char file_name[])
{
	size_t fw_ret;
	size_t max_nmemb = cal_max_nmemb(size);

	if(size == (size_t)0) return(NORMAL_RC);
	if(nmemb == (size_t)0) return(NORMAL_RC);

	if((fp == NULL)||(ptr == NULL)) return(ARG_ERROR_RC);

	while(nmemb > max_nmemb){
		fw_ret = fwrite(ptr, size, max_nmemb, fp);
		if(fw_ret != max_nmemb){
			log_printf(1, "\n%d*%d Bytes write failed.\n",
			           (int)size, (int)max_nmemb);
			if(file_name != NULL){
				RC_TRY( error_printf(10007, file_name) );
			}else{
				RC_TRY( error_printf(10006) );
			}
			return(WRITE_ERROR_RC);
		}
		nmemb -= max_nmemb;
		ptr = (void *)(((char *)ptr) + (size*max_nmemb));
	}

	fw_ret = fwrite(ptr, size, nmemb, fp);
	if(fw_ret != nmemb){
		log_printf(1, "\n%d*%d Bytes write failed.\n", (int)size, (int)nmemb);
		if(file_name != NULL){
			RC_TRY( error_printf(10007, file_name) );
		}else{
			RC_TRY( error_printf(10006) );
		}
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


/* fread() の RC 型ラッパー．*/
RC
rc_fread(void *ptr, size_t size, size_t nmemb, FILE *fp)
{
	RC_TRY( rc_fread_fn(ptr, size, nmemb, fp, NULL) );

	return(NORMAL_RC);
}


/* fread() の RC 型ラッパー．
 * エラー時にファイル名を表示する．
 */
RC
rc_fread_fn(void *ptr, size_t size, size_t nmemb, FILE *fp,
            const char file_name[])
{
	size_t fr_ret;
	size_t max_nmemb = cal_max_nmemb(size);

	if(size == (size_t)0) return(NORMAL_RC);
	if(nmemb == (size_t)0) return(NORMAL_RC);

	if((fp == NULL)||(ptr == NULL)) return(ARG_ERROR_RC);

	while(nmemb > max_nmemb){
/* fprintf(stderr, "r:%d %d", size, max_nmemb); */
		fr_ret = fread(ptr, size, max_nmemb, fp);
/* fprintf(stderr, " done.\n"); */
		if(fr_ret != max_nmemb){
			log_printf(1, "\n%d*%d Bytes read failed.\n",
			           (int)size, (int)max_nmemb);
			if(file_name != NULL){
				RC_TRY( error_printf(10009, file_name) );
			}else{
				RC_TRY( error_printf(10008) );
			}
			return(WRITE_ERROR_RC);
		}
		nmemb -= max_nmemb;
		ptr = (void *)(((char *)ptr) + (size*max_nmemb));
	}

/* fprintf(stderr, "r:%d %d", size, nmemb); */
	fr_ret = fread(ptr, size, nmemb, fp);
/* fprintf(stderr, " done.\n"); */
	if(fr_ret != nmemb){
		log_printf(1, "\n%d*%d Bytes read failed.\n", (int)size, (int)nmemb);
		if(file_name != NULL){
			RC_TRY( error_printf(10009, file_name) );
		}else{
			RC_TRY( error_printf(10008) );
		}
		return(WRITE_ERROR_RC);
	}

	return(NORMAL_RC);
}


