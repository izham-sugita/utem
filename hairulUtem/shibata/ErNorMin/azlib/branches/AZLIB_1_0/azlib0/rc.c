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

/* $Id: rc.c,v 1.5 2003/07/22 10:21:27 adachi Exp $ */


#include <stdio.h>
#include <stdlib.h>
#include "rc.h"

#define FUNC_SIZE 10

/* 引数を４つ持つ関数のポインタの配列 */
static void (*RcErrorPrintFunc[FUNC_SIZE])(const char *, int,
                                           const char *, RC *)
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


RC rc_malloc(size_t size, void **ptr)
{
	if(size <= 0){
		return(ARG_ERROR_RC);
	}

	(*ptr) = malloc(size);
	if(*ptr == NULL){
		return(ALLOC_ERROR_RC);
	}

	return(NORMAL_RC);
}


RC rc_fopen(const char *path, const char *mode, FILE **fp)
{
	*fp = fopen(path, mode);
	if(*fp == NULL){
		RC rc = OPEN_ERROR_RC;
		rc_error_print(__FILE__, __LINE__, path, &rc);
		return(rc);
	}

	return(NORMAL_RC);
}


RC rc_fclose(FILE *fp)
{
	int ret;

	ret = fclose(fp);
	if(ret != 0) return(CLOSE_ERROR_RC);

	return(NORMAL_RC);
}


void print_rc(FILE *stream, RC rc)
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
		fprintf(stream, "CONVERT_ERROR_RC");
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


/* エラー出力（デフォルト） */
void rc_error_print_default(const char *file_name, int line,
                            const char *rc_func, RC *rc)
{
	fprintf(stderr, "[%s : %d] %s\n", file_name, line, rc_func);
	fprintf(stderr, "    ==> ");
	print_rc(stderr, *rc);
}


/* エラー出力（何もしない） */
void rc_error_print_silent(const char *file_name, int line,
                           const char *rc_func, RC *rc)
{
	return;
}


/* エラー出力関数を func に変更 */
void rc_push_func(void (*func)(const char *, int, const char *, RC *))
{
	StackPtr++;
	if(StackPtr >= FUNC_SIZE) StackPtr = FUNC_SIZE - 1;

	if(func == NULL) func = rc_error_print_default;
	RcErrorPrintFunc[StackPtr] = func;
}


/* エラー出力関数を元に戻す */
void rc_pop_func(void)
{
	StackPtr--;
	if(StackPtr < 0) StackPtr = 0;
}


/* エラー出力関数ドライバ */
void rc_error_print(const char *file_name, int line,
                    const char *rc_func, RC *rc)
{
	(RcErrorPrintFunc[StackPtr])(file_name, line, rc_func, rc);
}


