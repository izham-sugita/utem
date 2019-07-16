/*********************************************************************
 * log_printf.c
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *  <Takuya HAYASHI>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: log_printf.c 1043 2016-02-28 23:58:30Z hayashi $ */

#include <stdio.h>
#include <stdarg.h>
#include <time.h>
#include "rc.h"
#include "time_utl.h"

#define MAX_ID_COUNT    (10)

typedef struct{
	int id;
	char *str;
} STRING_TABLE;

/* エラー，ワーニングメッセージ */
#include "messages.h"

/* ログファイルのポインタ */
static FILE *LogFile[] = {NULL, NULL, NULL, NULL, NULL};
/* ログレベル */
static int LogLevel = 3;

static int search_string_table(int table_size,
                               const STRING_TABLE table[], int id);

/* RC型の関数 rc_func を実行したときの返り値が NORMAL_RC 以外のときに呼び出さ
 * れる関数で，message.h で定義されたエラーメッセージおよびエラーが発生した部分
 * のファイル名と行番号をすべてのログ出力先に出力する．
 *
 * rc_push_func()でエラー出力関数として設定して使うため，直接呼び出すことはな
 * い．
 */
void
rc_error_print_log(const char *file_name, int line,
                   const char *rc_func, RC rc)
{
	int ii1;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int index = search_string_table(sizeof(ErrorStr)/sizeof(ErrorStr[0]),
	                                ErrorStr, 99999);
	static int flag = 0;

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			if( (flag == 0)&&(index >= 0) ){
				fprintf(LogFile[ii1], "\n| ERROR %5.5d: ", 99999);
				fprintf(LogFile[ii1], "%s\n", ErrorStr[index].str);
			}
			fprintf(LogFile[ii1], "\n[%s : %d] ", file_name, line);
			print_rc(LogFile[ii1], rc);
		}
	}
	flag = 1;
}

/* ログ出力レベルを level に設定する．
 */
RC
set_log_level(int level)
{
	LogLevel = level;
	return(NORMAL_RC);
}

/* 現在のログ出力レベルを返す．
 */
int
get_log_level(void)
{
	return(LogLevel);
}

/* i 番目のログ出力先を fp に設定する．
 */
RC
set_log_file(int i, FILE *fp)
{
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(i < 0) return(ARG_ERROR_RC);
	if(i >= log_file_size) return(ARG_ERROR_RC);

	LogFile[i] = fp;

	return(NORMAL_RC);
}

/* i 番目のログ出力先のファイルポインタを返す．
 */
FILE *
get_log_file(int i)
{
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(i < 0) return(NULL);
	if(i >= log_file_size) return(NULL);

	return(LogFile[i]);
}

/* ログ出力レベル level で，書式 format, ... に従った文字列を，タイムスタンプ
 * 付きですべてのログファイルに出力する．
 */
RC
tlog_printf(int level, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	char buf[128];


	if(level > LogLevel) return(NORMAL_RC);

	RC_TRY( time_string_stamp(buf, sizeof(buf)) );

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			if(fprintf(LogFile[ii1], "%s : ", buf) < 0) return(WRITE_ERROR_RC);
			/* fflush(LogFile[ii1]); */
			va_start(ap, format);
			if(vfprintf(LogFile[ii1], format, ap) < 0) return(WRITE_ERROR_RC);
			/* fflush(LogFile[ii1]); */
			va_end(ap);
		}
	}

	return(NORMAL_RC);
}


/* ログ出力レベル level で，書式 format, ... に従った文字列をすべてのログファイ
 * ルに出力する．
 */
RC
log_printf(int level, const char *format, ...)
{
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, format);
			if(vfprintf(LogFile[ii1], format, ap) < 0) return(WRITE_ERROR_RC);
			/* fflush(LogFile[ii1]); */
			va_end(ap);
		}
	}

	return(NORMAL_RC);
}

#define LINE_LENGTH (80)
/* 全体で半角 LINE_LENGTH 字になるように横線を付けたlog_printf．
 * 横線に使う文字をline_charで設定する．
 */
RC
log_printf_line(int level, char line_char, const char *format, ...)
{
	int ii1;
	int size_mes = 0;
	int size_line;
	char sub[LINE_LENGTH+2];
	char message[LINE_LENGTH+3];
	va_list ap; 
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(level > LogLevel) return(NORMAL_RC);

	va_start(ap, format);
	if(vsprintf(sub, format, ap) < 0) return(WRITE_ERROR_RC);
	/* format, ... の長さを計算する． */
	while(sub[size_mes] != '\0'){
		size_mes++;
	}   
	/* 横線を付ける前に LINE_MESSAGE 字を超える場合は，エラーとする． */
	if(size_mes > LINE_LENGTH){
		fprintf(stderr, "Make the message shorter.\n");
		return(ARG_ERROR_RC);
	}   
	size_line = LINE_LENGTH - 2 - size_mes;

	/* メッセージの作成 */
	message[0] = '\n';
	for(ii1=1; ii1<size_line/2-1; ii1++){
		message[ii1] = line_char;
	}
	message[size_line/2-1] = ' ';
	for(ii1=0; ii1<size_mes; ii1++){
		message[size_line/2 + ii1] = sub[ii1];
	}
	message[size_line/2 + size_mes] = ' ';
	for(ii1=size_line/2 + size_mes + 1; ii1<80; ii1++){
		message[ii1] = line_char;
	}
	message[LINE_LENGTH] = '\n';
	message[LINE_LENGTH+1] = '\0';

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			if(fprintf(LogFile[ii1], message) < 0) return(WRITE_ERROR_RC);
		}
	}
	va_end(ap);

	return(NORMAL_RC);
}


/* ログ出力レベル level で，書式 format に従った文字列を，以前出力した分を消去し
 * てすべてのログファイルに出力する．パーセンテージなどの進捗状況を同じ行に出力
 * したいときに使う．
 * 
 * コンソール（標準出力および標準エラー出力）にしか出力されない．
 */
RC
progress_printf(int level, const char *format, ...)
{
	int ii1;
	int bs;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(level > LogLevel) return(NORMAL_RC);

	for(ii1=0; ii1<log_file_size; ii1++){
		/* don't output files. display only!!  */
		if((LogFile[ii1] == stderr)||(LogFile[ii1] == stdout)){
			bs = 0;
			va_start(ap, format);
			if( (bs = vfprintf(LogFile[ii1], format, ap)) < 0 )
				return(WRITE_ERROR_RC);
			va_end(ap);
			for(; bs; bs--) if(fputc('\b', LogFile[ii1]) == EOF)
				return(WRITE_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}


/* 進捗状況のパーセンテージをプログレスバー付きで表示する．
 * パーセンテージの値は 0.0以上，1.0以下
 * 例: [======>    ] 50% 123 hoge
 */
#define PROGRESS_BAR_SIZE (30)
RC
progress_bar_printf(int level, double percent, const char *format, ...)
{
	int ii1, ii2;
	int bs;
	va_list ap;
	char buf[PROGRESS_BAR_SIZE+1+4];
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	if(level > LogLevel) return(NORMAL_RC);

	if(percent < 0.0) percent = 0.0;
	if(percent > 1.0) percent = 1.0;

	for(ii1=0; ii1<log_file_size; ii1++){
		/* don't output files. display only!!  */
		if((LogFile[ii1] == stderr)||(LogFile[ii1] == stdout)){
			/* [======>    ] 50% */
			buf[0] = '[';
			for(ii2=1; ii2<(percent*PROGRESS_BAR_SIZE-1); ii2++) buf[ii2] = '=';
			if(ii2 < PROGRESS_BAR_SIZE-1) buf[ii2++] = '>';
			for(; ii2<PROGRESS_BAR_SIZE-1; ii2++) buf[ii2] = ' ';
			buf[PROGRESS_BAR_SIZE-1] = ']';
			if( sprintf(&buf[++ii2], "%3d%%", (int)(percent*100)) != 4 )
				return(WRITE_ERROR_RC);
			if( fputs(buf, LogFile[ii1]) == EOF) return(WRITE_ERROR_RC);

			/* プログレスバー右端の表示 "[======>    ] 50%(ここ)" */
			bs = 0;
			if(format != NULL){
				va_start(ap, format);
				if( (bs = vfprintf(LogFile[ii1], format, ap)) < 0 )
				   	return(WRITE_ERROR_RC);
				va_end(ap);
			}
 
			/* 後退 */
			bs += (PROGRESS_BAR_SIZE+4);
			for(; bs; bs--)
				if(fputc('\b', LogFile[ii1]) == EOF) return(WRITE_ERROR_RC);
			fflush(LogFile[ii1]);
		}
	}

	return(NORMAL_RC);
}


/* すべてのログファイルをフラッシュする．
 */
RC
log_flush(void)
{
	int ii1;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);

	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL) fflush(LogFile[ii1]);
	}

	return(NORMAL_RC);
}

/* エラーメッセージの配列の中から，エラーIDで指定したメッセージのインデックスを返す．
 */
static int
search_string_table(int table_size, const STRING_TABLE table[], int id)
{
	int ii1;

	for(ii1=0; ii1<table_size; ii1++){
		if(table[ii1].id == id) return(ii1);
	}

	return(-1);
}


/* messages.h のエラーID error_id の書式に従ったエラーメッセージを出力する．
 * MAX_ID_COUNT で指定されている回数以降は出力しない．
 */
RC
error_printf(int error_id, ...)
{
	static int id_count[][2] = {{-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}};
	static int id_count_index = 0;
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int index = search_string_table(sizeof(ErrorStr)/sizeof(ErrorStr[0]),
	                                ErrorStr, error_id);
	RC_NEG_CHK( index );

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
			fprintf(LogFile[ii1], "\n| ERROR %5.5d: ", error_id);
			if(vfprintf(LogFile[ii1], ErrorStr[index].str, ap) < 0){
				return(WRITE_ERROR_RC);
			}
			fflush(LogFile[ii1]);
			va_end(ap);
		}
	}

	rc_push_func(rc_error_print_silent);

	return(NORMAL_RC);
}


/* messages.h のエラーID error_id の書式に従ったワーニングメッセージを出力する．
 * MAX_ID_COUNT で指定されている回数以降は出力しない．
 */
RC
warning_printf(int warning_id, ...)
{
	static int id_count[][2] = {{-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}, {-1, 0}};
	static int id_count_index = 0;
	int ii1;
	va_list ap;
	int log_file_size = sizeof(LogFile)/sizeof(LogFile[0]);
	int index = search_string_table(sizeof(WarningStr)/sizeof(WarningStr[0]),
	                                WarningStr, warning_id);
	RC_NEG_CHK( index );

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

	RC_NEG_CHK( index );
	for(ii1=0; ii1<log_file_size; ii1++){
		if(LogFile[ii1] != NULL){
			va_start(ap, warning_id);
			fprintf(LogFile[ii1], "\n| WARNING %5.5d: ", warning_id);
/*			{
				wchar_t wbuf[1024];
				mbstowcs(wbuf, WarningStr[index].str,
				         strlen(WarningStr[index].str) + 1);
				if(vfwprintf(LogFile[ii1], wbuf, ap) < 0){
					return(WRITE_ERROR_RC);
				}
			}*/
			if(vfprintf(LogFile[ii1], WarningStr[index].str, ap) < 0){
				return(WRITE_ERROR_RC);
			}
			fflush(LogFile[ii1]);
			va_end(ap);
		}
	}

	return(NORMAL_RC);
}

