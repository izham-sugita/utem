/*********************************************************************
 * log_printf.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: log_printf.c 645 2005-12-09 02:23:38Z sasaoka $ */

#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include "log_printf.h"
#include "rc.h"
#include "time_utl.h"

#define MAX_ID_COUNT    (10)

typedef struct{
	int id;
	char *str;
} STRING_TABLE;

/* エラー，ワーニングメッセージ */
#include "messages.h"

static FILE *LogFile[] = {NULL, NULL, NULL, NULL, NULL};
static int LogLevel = 3;

static int search_string_table(int table_size,
                               const STRING_TABLE table[], int id);
static void print_filename(FILE *fp, const char *file_name);
static void rc_error_print_log_sub(FILE *fp, const char *file_name, int line,
                                   const char *rc_func, RC rc);
static RC log_printf_default(FILE *fp, const char *format, va_list ap);
static RC tlog_printf_default(FILE *fp, const char *format, va_list ap);
static RC progress_printf_default(FILE *fp, const char *format, va_list ap);
static RC progress_bar_printf_default(FILE *fp, double percent,
                                      const char *format, va_list ap);
static RC error_printf_default(FILE *fp, int error_id, va_list ap);
static RC warning_printf_default(FILE *fp, int warning_id, va_list ap);

static RC (*LogPrintfFunc)
          (FILE *, const char *, va_list) = log_printf_default;
static RC (*TLogPrintfFunc)
          (FILE *, const char *, va_list) = tlog_printf_default;
static RC (*ProgressPrintfFunc)
          (FILE *, const char *, va_list) = progress_printf_default;
static RC (*ProgressBarPrintfFunc)
          (FILE *, double, const char *,va_list) = progress_bar_printf_default;
static RC (*ErrorPrintfFunc)
          (FILE *, int, va_list) = error_printf_default;
static RC (*WarningPrintfFunc)
          (FILE *, int, va_list) = warning_printf_default;


void
rc_error_print_log (const char *file_name, int line,
                    const char *rc_func, RC rc)
{
	int ii1;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			rc_error_print_log_sub(LogFile[ii1], file_name, line, rc_func, rc);
			count++;
		}
	}
	if(count == 0){
		rc_error_print_log_sub(stderr, file_name, line, rc_func, rc);
	}
}


static void
rc_error_print_log_sub (FILE *fp, const char *file_name, int line,
                        const char *rc_func, RC rc)
{
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	static int index = -1;
	static int flag = 0;

	if(flag == 0){
		index = search_string_table(sizeof(ErrorStr)/sizeof(ErrorStr[0]),
	                                ErrorStr, 99999);
	}
	if( (flag < log_file_size)&&(index >= 0) ){
		fprintf(fp, "\n| ERROR %5.5d: ", 99999);
		fprintf(fp, "%s\n", ErrorStr[index].str);
	}
	fprintf(fp, "\n[%s : %d] ", file_name, line);
	print_rc(fp, rc);
	print_filename(fp, file_name);
	fprintf(fp, ": %4.4X : %2.2X\n", line, (int)rc);
	fflush(fp);
	flag++;
}


static void
print_filename (FILE *fp, const char *file_name)
{
	int ii1;

	for(ii1=0; file_name[ii1]; ii1++){
		fprintf(fp, "%2.2X ", file_name[ii1]);
	}
	for(; ii1<16; ii1++){
		fprintf(fp, "%2.2X ", 0);
	}
}


RC
set_log_level (int level)
{
	LogLevel = level;
	return(NORMAL_RC);
}


int
get_log_level (void)
{
	return(LogLevel);
}


RC
set_log_file (int i, FILE *fp)
{
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(i < 0) return(ARG_ERROR_RC);
	if(i >= log_file_size) return(ARG_ERROR_RC);

	LogFile[i] = fp;

	return(NORMAL_RC);
}


RC
log_printf (int level, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, format);
			RC_TRY( LogPrintfFunc(LogFile[ii1], format, ap) );
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, format);
		RC_TRY( LogPrintfFunc(stderr, format, ap) );
		va_end(ap);
	}

	return(NORMAL_RC);
}


static RC
log_printf_default (FILE *fp, const char *format, va_list ap)
{
	if(fp != NULL){
		if(vfprintf(fp, format, ap) < 0) return(WRITE_ERROR_RC);
		fflush(fp);
	}
	return(NORMAL_RC);
}

void
set_log_printf_func (RC (*func)(FILE *, const char *, va_list))
{
	if(func != NULL) LogPrintfFunc = func;
}

void
reset_log_printf_func (void)
{
	LogPrintfFunc = log_printf_default;
}


RC
tlog_printf (int level, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, format);
			RC_TRY( TLogPrintfFunc(LogFile[ii1], format, ap) );
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, format);
		RC_TRY( TLogPrintfFunc(stderr, format, ap) );
		va_end(ap);
	}

	return(NORMAL_RC);
}


static RC
tlog_printf_default (FILE *fp, const char *format, va_list ap)
{
	char buf[128];

	RC_TRY( time_string_stamp(buf, sizeof(buf)) );

	if(fp != NULL){
		if(fprintf(fp, "%s : ", buf) < 0) return(WRITE_ERROR_RC);
		if(vfprintf(fp, format, ap) < 0) return(WRITE_ERROR_RC);
		fflush(fp);
	}

	return(NORMAL_RC);
}


void
set_tlog_printf_func (RC (*func)(FILE *, const char *, va_list))
{
	if(func != NULL) TLogPrintfFunc = func;
}

void
reset_tlog_printf_func (void)
{
	TLogPrintfFunc = tlog_printf_default;
}


RC
progress_printf (int level, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		/* don't output files  */
		if((LogFile[ii1] == stderr)||(LogFile[ii1] == stdout)){
			va_start(ap, format);
			RC_TRY( ProgressPrintfFunc(LogFile[ii1], format, ap) );
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, format);
		RC_TRY( ProgressPrintfFunc(stderr, format, ap) );
		va_end(ap);
	}

	return(NORMAL_RC);
}


static RC
progress_printf_default (FILE *fp, const char *format, va_list ap)
{
	int bs;
	if(fp != NULL){
		if( (bs = vfprintf(fp, format, ap)) < 0 ) return(WRITE_ERROR_RC);
		for(; bs; bs--) if(fputc('\b', fp) == EOF) return(WRITE_ERROR_RC);
		fflush(fp);
	}
	return(NORMAL_RC);
}

void
set_progress_printf_func (RC (*func)(FILE *, const char *, va_list))
{
	if(func != NULL) ProgressPrintfFunc = func;
}

void
reset_progress_printf_func (void)
{
	ProgressPrintfFunc = progress_printf_default;
}


RC
progress_bar_printf (int level, double percent, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		/* don't output files  */
		if((LogFile[ii1] == stderr)||(LogFile[ii1] == stdout)){
			va_start(ap, format);
			RC_TRY( ProgressBarPrintfFunc(LogFile[ii1], percent, format, ap) );
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, format);
		RC_TRY( ProgressBarPrintfFunc(stderr, percent, format, ap) );
		va_end(ap);
	}

	return(NORMAL_RC);
}


#define PROGRESS_BAR_SIZE (30)
static RC
progress_bar_printf_default (FILE *fp, double percent,
                             const char *format, va_list ap)
{
	int ii1, bs = 0;
	char buf[PROGRESS_BAR_SIZE+1+4];

	if(percent < 0.0) percent = 0.0;
	if(percent > 1.0) percent = 1.0;

	/* [======>    ] 50% */
	buf[0] = '[';
	for(ii1=1; ii1<(percent*PROGRESS_BAR_SIZE-1); ii1++) buf[ii1] = '=';
	if(ii1 < PROGRESS_BAR_SIZE-1) buf[ii1++] = '>';
	for(; ii1<PROGRESS_BAR_SIZE-1; ii1++) buf[ii1] = ' ';
	buf[PROGRESS_BAR_SIZE-1] = ']';
	if( sprintf(&buf[++ii1], "%3d%%", (int)(percent*100)) != 4 )
		return(WRITE_ERROR_RC);
	if( fputs(buf, fp) == EOF) return(WRITE_ERROR_RC);

	/* プログレスバー右の表示 "[======>    ] 50%(ここ)" */
	if(format != NULL){
		if( (bs = vfprintf(fp, format, ap)) < 0 ) return(WRITE_ERROR_RC);
	}

	/* 後退 */
	bs += (PROGRESS_BAR_SIZE+4);
	for(; bs; bs--) if(fputc('\b', fp) == EOF) return(WRITE_ERROR_RC);
	fflush(fp);

	return(NORMAL_RC);
}

void
set_progress_bar_printf_func (RC (*func)(FILE *, double percent,
                              const char *, va_list))
{
	if(func != NULL) ProgressBarPrintfFunc = func;
}

void
reset_progress_bar_printf_func (void)
{
	ProgressBarPrintfFunc = progress_bar_printf_default;
}


/* 全てのログファイルをフラッシュ */
RC
log_flush (void)
{
	int ii1;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL) fflush(LogFile[ii1]);
	}

	return(NORMAL_RC);
}


static int
search_string_table(int table_size, const STRING_TABLE table[], int id)
{
	int ii1;

	for(ii1=0; ii1<table_size; ii1++){
		if(table[ii1].id == id) return(ii1);
	}

	return(-1);
}


RC
error_printf (int error_id, ...)
{
	static int id_count[][2] = {{-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}};
	static int id_count_index = 0;
	int ii1;
	va_list ap;
	int count = 0;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(LogLevel < 0) return(NORMAL_RC);
	for(ii1=0; ii1<(sizeof(id_count)/sizeof(id_count[0])); ii1++){
		if(id_count[ii1][0] == error_id){
			if(id_count[ii1][1] >= MAX_ID_COUNT){
				if(id_count[ii1][1] == MAX_ID_COUNT){
					RC_TRY( warning_printf(90001, error_id) );
					id_count[ii1][1]++;
				}
				return(NORMAL_RC);
			}else{
				id_count[ii1][1]++;
			}
			break;
		}
	}
	if(ii1 >= (sizeof(id_count)/sizeof(id_count[0]))){
		id_count[id_count_index][0] = error_id;
		id_count[id_count_index][1] = 1;
		id_count_index++;
		if(id_count_index >= (sizeof(id_count)/sizeof(id_count[0]))){
			id_count_index = 0;
		}
	}

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, error_id);
			RC_TRY( ErrorPrintfFunc(LogFile[ii1], error_id, ap) );
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, error_id);
		RC_TRY( ErrorPrintfFunc(stderr, error_id, ap) );
		va_end(ap);
	}

	rc_push_func(rc_error_print_silent);

	return(NORMAL_RC);
}


static RC
error_printf_default (FILE *fp, int error_id, va_list ap)
{
	int index = -1;
	static int flag = 0;

	if(flag == 0){
		index = search_string_table(sizeof(ErrorStr)/sizeof(ErrorStr[0]),
		                            ErrorStr, error_id);
		RC_NEG_CHK( index );
		flag = 1;
	}
	fprintf(fp, "\n| ERROR %5.5d: ", error_id);
	if(vfprintf(fp, ErrorStr[index].str, ap) < 0) return(WRITE_ERROR_RC);
	fflush(fp);
	return(NORMAL_RC);
}

void
set_error_printf_func (RC (*func)(FILE *, int , va_list))
{
	if(func != NULL)  ErrorPrintfFunc = func;
}

void
reset_error_printf_func (void)
{
	ErrorPrintfFunc = error_printf_default;
}


RC
warning_printf (int warning_id, ...)
{
	static int id_count[][2] = {{-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}};
	static int id_count_index = 0;
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int count = 0;

	if(LogLevel < 0) return(NORMAL_RC);
	for(ii1=0; ii1<(sizeof(id_count)/sizeof(id_count[0])); ii1++){
		if(id_count[ii1][0] == warning_id){
			if(id_count[ii1][1] >= MAX_ID_COUNT){
				if(id_count[ii1][1] == MAX_ID_COUNT){
					RC_TRY( warning_printf(90002, warning_id) );
					id_count[ii1][1]++;
				}
				return(NORMAL_RC);
			}else{
				id_count[ii1][1]++;
			}
			break;
		}
	}
	if(ii1 >= (sizeof(id_count)/sizeof(id_count[0]))){
		id_count[id_count_index][0] = warning_id;
		id_count[id_count_index][1] = 1;
		id_count_index++;
		if(id_count_index >= (sizeof(id_count)/sizeof(id_count[0]))){
			id_count_index = 0;
		}
	}

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, warning_id);
			RC_TRY( WarningPrintfFunc(LogFile[ii1], warning_id, ap) );
			fflush(LogFile[ii1]);
			va_end(ap);
			count++;
		}
	}
	if(count == 0){
		va_start(ap, warning_id);
		RC_TRY( WarningPrintfFunc(stderr, warning_id, ap) );
		va_end(ap);
	}

	return(NORMAL_RC);
}


static RC
warning_printf_default (FILE *fp, int warning_id, va_list ap)
{
	int index = -1;
	static int flag = 0;

	if(flag == 0){
		index = search_string_table(sizeof(WarningStr)/sizeof(WarningStr[0]),
		                            WarningStr, warning_id);
		RC_NEG_CHK( index );
		flag = 1;
	}

	fprintf(fp, "\n| WARNING %5.5d: ", warning_id);
	if(vfprintf(fp, WarningStr[index].str, ap) < 0) return(WRITE_ERROR_RC);
	fflush(fp);

	return(NORMAL_RC);
}

void
set_warning_printf_func (RC (*func)(FILE *, int , va_list))
{
	if(func != NULL)  WarningPrintfFunc = func;
}

void
reset_warning_printf_func (void)
{
	WarningPrintfFunc = warning_printf_default;
}



