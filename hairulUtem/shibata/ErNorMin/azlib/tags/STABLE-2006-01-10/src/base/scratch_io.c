/*********************************************************************
 * scratch_io.c
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

/* $Id: scratch_io.c 655 2005-12-09 07:42:03Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#ifndef WIN32
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#endif  /* WIN32 */
#include "scratch_io.h"
#include "wrapper_64.h"
#include "log_printf.h"
#include "macros.h"


#define SF_ARRAY_SIZE   (10)      /* 使用可能なスクラッチファイルの最大値 */
static SCRATCH_FILE SFArray[SF_ARRAY_SIZE];
static int InitFlag = 0;

static void delete_all_scratch_file(void);


SCRATCH_FILE *
open_scratch_file (const char scratch_dir[])
{
	int ii1;
	SCRATCH_FILE *sf;

	if(InitFlag == 0){
		for(ii1=0; ii1<SF_ARRAY_SIZE; ii1++){
			SFArray[ii1].fp = NULL;
			SFArray[ii1].file_name = NULL;
		}
		InitFlag = 1;
		atexit(delete_all_scratch_file);
	}

	sf = NULL;
	for(ii1=0; ii1<SF_ARRAY_SIZE; ii1++){
		if(SFArray[ii1].fp == NULL){
			sf = &(SFArray[ii1]);
			break;
		}
	}
	if(sf == NULL) return(NULL);

#ifdef WIN32
	if((sf->file_name = _tempnam(NULL, "__scratch_work_file__.")) == NULL){
		return(NULL);
	}
	if((sf->fp = fopen(sf->file_name, "w+bc")) == NULL) return(NULL);
#else  /* WIN32 */
	{
		int fd;
		mode_t old_mode;
		/* XXXXXXを変更しないこと．See man mkstemp. */
		char tmpfname_base[] = "__scratch_work_file__XXXXXX";
		int str_len;

		str_len = strlen(tmpfname_base) + 1;
		if(scratch_dir) str_len += strlen(scratch_dir) + 1;

		if((sf->file_name = malloc(str_len)) == NULL) return(NULL);
		if(scratch_dir){
			sprintf(sf->file_name, "%s/%s",
			        scratch_dir, tmpfname_base);
		}else{
			sprintf(sf->file_name, "%s", tmpfname_base);
		}

		old_mode = umask(077); /* 600(所有者のみ読書可)に設定 */
		fd = mkstemp(sf->file_name);
		umask(old_mode); /* 以前のパーミッションモードに戻す */
		if(fd == -1) return(NULL);
		if(!(sf->fp = fdopen(fd, "w+bc"))) return(NULL);
		unlink(sf->file_name); /* ファイルを見えなくする */
	}
#endif /* WIN32 */
	tlog_printf(5, "Scratch file [%s] is open.\n", sf->file_name);

	for(ii1=0; ii1<SCRATCH_DATA_SIZE; ii1++){
		sf->position[ii1] = (WRAP64_INT)0;
		sf->use_flag[ii1] = 0;
	}
	sf->stack_ptr = 0;

	return(sf);
}


/* atexit() 用 */
static void
delete_all_scratch_file (void)
{
	int ii1;

	for(ii1=0; ii1<SF_ARRAY_SIZE; ii1++){
		if(SFArray[ii1].fp != NULL){
			fclose(SFArray[ii1].fp);
			SFArray[ii1].fp = NULL;
			if(SFArray[ii1].file_name != NULL){
				remove(SFArray[ii1].file_name);
				SFArray[ii1].file_name = NULL;
			}
		}
	}
}


RC
close_scratch_file (SCRATCH_FILE *sf)
{
	int ii1;

	if(sf->fp == NULL) return(ARG_ERROR_RC);
	if(fclose(sf->fp) != 0) return(CLOSE_ERROR_RC);
	sf->fp = NULL;
#ifdef WIN32
	if(remove(sf->file_name) != 0) return(UNKNOWN_ERROR_RC);
#endif /* WIN32 */
	RC_TRY( tlog_printf(5, "Scratch file [%s] is closed.\n", sf->file_name) );
	free(sf->file_name);
	sf->file_name = NULL;

	for(ii1=0; ii1<SCRATCH_DATA_SIZE; ii1++){
		sf->position[ii1] = (WRAP64_INT)0;
		sf->use_flag[ii1] = 0;
	}
	sf->stack_ptr = 0;

	return(NORMAL_RC);
}


/* スクラッチファイルへの書き込みは，この関数の戻り値を使って */
/* シーケンシャルに行うこと */
FILE *
add_scratch_file_start (SCRATCH_FILE *sf)
{
	if(sf->fp == NULL) return(NULL);
	if(sf->stack_ptr >= SCRATCH_DATA_SIZE) return(NULL);
/*	fflush(sf->fp); */
	if(wrap64_fseek(sf->fp, sf->position[sf->stack_ptr], SEEK_SET)
	                                                  != NORMAL_RC){
		return(NULL);
	}
	sf->use_flag[sf->stack_ptr] = 2;   /* 書き込み中 */

	return(sf->fp);
}


/* スクラッチファイルへの書き込みが終了したら，この関数を呼び出す */
/* 戻り値は，読み込みと削除の際に使用する */
int
add_scratch_file_end (SCRATCH_FILE *sf)
{
	int ret;

	if(sf->fp == NULL) return(-1);
	if(sf->stack_ptr >= SCRATCH_DATA_SIZE) return(-2);
	if(sf->use_flag[sf->stack_ptr] != 2) return(-3);

	sf->use_flag[sf->stack_ptr] = 1;   /* 使用中 */
/*	fflush(sf->fp); */
	ret = sf->stack_ptr;
	(sf->stack_ptr)++;
	sf->position[sf->stack_ptr] = wrap64_ftell(sf->fp);

	return(ret);
}


/* スクラッチファイルの読み込みは，この関数の戻り値を使用する */
/* handle : add_scratch_file_end() の戻り値 */
FILE *
read_scratch_file (SCRATCH_FILE *sf, int handle)
{
	if(sf->fp == NULL) return(NULL);
	if(handle >= sf->stack_ptr) return(NULL);
	if(sf->use_flag[handle] != 1) return(NULL);
/*	fflush(sf->fp); */
	if(wrap64_fseek(sf->fp, sf->position[handle], SEEK_SET) != NORMAL_RC){
		return(NULL);
	}

	return(sf->fp);
}


/* 既存データの上書き */
/* handle : add_scratch_file_end() の戻り値 */
FILE *
write_scratch_file (SCRATCH_FILE *sf, int handle)
{
	return(read_scratch_file(sf, handle)); /* 実は読み込みと同じ処理 */
}


/* add_scratch_file_end() の戻り値 handle のデータを削除する */
RC
delete_scratch_file (SCRATCH_FILE *sf, int handle)
{
	int ii1;

	if(sf->fp == NULL) return(ARG_ERROR_RC);
	if(sf->stack_ptr <= 0) return(ARG_ERROR_RC);
	if(handle >= sf->stack_ptr) return(ARG_ERROR_RC);
	if(sf->use_flag[handle] != 1) return(ARG_ERROR_RC);

	sf->use_flag[handle] = 0;

	for(ii1=(sf->stack_ptr-1); ii1>=0; ii1--){
		if(sf->use_flag[ii1] != 0) break;
		sf->position[sf->stack_ptr] = 0L;
		(sf->stack_ptr)--;
	}

	return(NORMAL_RC);
}

RC
print_scratch_file_info (FILE *fp, SCRATCH_FILE *sf)
{
	int ii1;

	if(sf == NULL) return(ARG_ERROR_RC);

	fprintf(fp, "fp = %p\n", sf->fp);
	fprintf(fp, "file_name = %s\n", sf->file_name);
	for(ii1=0; ii1<SCRATCH_DATA_SIZE; ii1++){
#ifdef WIN32
		fprintf(fp, "position[%3d] = %I64d use_flag = %s\n", ii1,
		        sf->position[ii1], (sf->use_flag[ii1]==1?"Yes":"No"));
#else
 #ifdef _LARGEFILE_SOURCE
		fprintf(fp, "position[%3d] = %lld use_flag = %s\n", ii1,
		        (long long)sf->position[ii1],
		        (sf->use_flag[ii1]==1?"Yes":"No"));
 #else
		fprintf(fp, "position[%3d] = %ld use_flag = %s\n", ii1,
		        (long)sf->position[ii1], (sf->use_flag[ii1]==1?"Yes":"No"));
 #endif
#endif
	}
	fprintf(fp, "stack_ptr = %d\n", sf->stack_ptr);

	return(NORMAL_RC);
}

