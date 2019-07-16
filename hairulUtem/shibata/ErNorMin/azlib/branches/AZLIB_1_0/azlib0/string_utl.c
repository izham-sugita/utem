/*********************************************************************
 * string_utl.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: string_utl.c,v 1.9 2004/02/05 07:15:02 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <time.h>
#include "rc.h"
#include "string_utl.h"


RC logprintf(char *f_name, char *format, ...)
{
	FILE *fp;
	va_list ap;
	struct tm *tm_ptr;
	time_t tt = time(NULL);
	char buf[128];

	tm_ptr = localtime(&tt);
	strftime(buf, sizeof(buf), "%d/%m/%Y", tm_ptr);
	buf[sizeof(buf)-1] = (char)'\0';

	RC_NULL_CHK( fp = fopen(f_name, "a") );
	va_start(ap, format);
	if(fprintf(fp, "%s : ", buf) < 0) return(WRITE_ERROR_RC);
	if(vfprintf(fp, format, ap) < 0) return(WRITE_ERROR_RC);

	va_end(ap);
	if(fclose(fp) != 0) return(CLOSE_ERROR_RC);

	return(NORMAL_RC);
}


RC init_string_array(STRING_ARRAY *stra)
{
	stra->size = 0;
	stra->alloc_size = 0;
	stra->str = NULL;

	return(NORMAL_RC);
}


RC free_string_array(STRING_ARRAY *stra)
{
	if(stra->alloc_size > 0){
		if(stra->str == NULL) return(ARG_ERROR_RC);

		free(stra->str);
	}
	stra->size = 0;
	stra->alloc_size = 0;
	stra->str = NULL;

	return(NORMAL_RC);
}


RC add_string_array(STRING_ARRAY *stra, const char *string)
{
	int length;
	
	if(string == NULL) return(ARG_ERROR_RC);
	if(stra->size > stra->alloc_size) return(ARG_ERROR_RC);
	if(stra->size == stra->alloc_size){
		(stra->alloc_size) += 64;
		stra->str = realloc(stra->str, (stra->alloc_size) * sizeof(char *));
		RC_NULL_CHK( stra->str );
	}
	RC_NEG_ZERO_CHK( length = strlen(string) + 1 );
	RC_NULL_CHK( stra->str[stra->size] = malloc(length * sizeof(char)) );
	strncpy(stra->str[stra->size], string, length);
	stra->str[stra->size][length-1] = (char)'\0';  /* 念の為 */
	(stra->size)++;

	return(NORMAL_RC);
}


void print_string_array(FILE *fp, STRING_ARRAY stra)
{
	int ii1;

	fprintf(fp, "size: %5d\n", stra.size);
	fprintf(fp, "alloc_size: %5d\n", stra.alloc_size);
	for(ii1=0;ii1<stra.size;ii1++){
		fprintf(fp, "string[%5d]: \"%s\"\n", ii1, stra.str[ii1]);
	}
}


RC split_string_delim(STRING_ARRAY *stra, const char *string, const char *delim)
{
	int ii1, ii2, ii3;
	int str_size, delim_flag;
	char *s_tmp;

	if(string == NULL) return(ARG_ERROR_RC);
	if(delim == NULL) return(ARG_ERROR_RC);
	str_size = (strlen(string)+1);
	RC_TRY( rc_malloc(str_size, (void *)&s_tmp) );
    for(ii1=0,ii2=0; ii1<str_size; ii1++){
		delim_flag = 0;
    	for(ii3=0; delim[ii3]; ii3++){
        	if(string[ii1] == delim[ii3]) delim_flag = 1;
		}
		if( (delim_flag == 1) || (string[ii1] == (char)'\0') ){
			if(ii2 != 0){
				s_tmp[ii2] = (char)'\0';
				add_string_array(stra, s_tmp);
				ii2 = 0;
			}
		}else{
			s_tmp[ii2] = string[ii1];
			ii2++;
		}
	}
	free(s_tmp);

	return(NORMAL_RC);
}


RC strip_string_char(char *string, const char *strip)
{
	char *p_tmp;
	int ii1, strip_flag;

	if(string == NULL) return(ARG_ERROR_RC);
	if(strip == NULL) return(ARG_ERROR_RC);
	for(p_tmp=string; *string; string++){
		strip_flag = 0;
		for(ii1=0; strip[ii1]; ii1++){
			if(*string == strip[ii1]) strip_flag = 1;
		}
		if(strip_flag != 1){
			*p_tmp++ = *string; 
		}
	}
	*p_tmp = '\0';

	return(NORMAL_RC);
}


RC print_backspace(int num)
{
	if(num < 0) return(NORMAL_RC);
	for(; num; num--) if(fputc('\b', stderr) == EOF) return(WRITE_ERROR_RC);
	return(NORMAL_RC);
}


RC print_progress(const char *format, ...)
{
	va_list ap;
	int bs;

	va_start(ap, format);
	if( (bs = vfprintf(stderr, format, ap)) < 0 ) return(WRITE_ERROR_RC);
	for(; bs; bs--) if(fputc('\b', stderr) == EOF) return(WRITE_ERROR_RC);
	va_end(ap);

	return(NORMAL_RC);
}


#define PROGRESS_BAR_SIZE (20)
RC print_progress_bar(double percent, const char *format, ...)
{
	int ii1, bs = 0;
	va_list ap;
	char buf[PROGRESS_BAR_SIZE+1];

	if(percent < 0.0) percent = 0.0;
	if(percent > 1.0) percent = 1.0;

	/* [======>50%    ] */
	buf[0] = '[';
	for(ii1=1; ii1<(percent*PROGRESS_BAR_SIZE-1); ii1++) buf[ii1] = '=';
	buf[ii1++] = '>';
	if(ii1 < PROGRESS_BAR_SIZE-3){
		if( sprintf(&buf[ii1], "%2d%%", (int)(percent*100)) < 0)
			return(WRITE_ERROR_RC);
		ii1+=3;
	}
	for(; ii1<PROGRESS_BAR_SIZE-1; ii1++) buf[ii1] = ' ';
	buf[PROGRESS_BAR_SIZE-1] = ']';
	buf[PROGRESS_BAR_SIZE] = '\0';
	if( fputs(buf, stderr) == EOF) return(WRITE_ERROR_RC);

	/* プログレスバー右の表示 */
	if(format != NULL){
		va_start(ap, format);
		if( (bs = vfprintf(stderr, format, ap)) < 0 ) return(WRITE_ERROR_RC);
		va_end(ap);
	}

	/* 後退 */
	bs += PROGRESS_BAR_SIZE;
	for(; bs; bs--) if(fputc('\b', stderr) == EOF) return(WRITE_ERROR_RC);

	return(NORMAL_RC);
}

