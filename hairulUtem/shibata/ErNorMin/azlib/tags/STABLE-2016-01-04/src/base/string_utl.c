/*********************************************************************
 * string_utl.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA> 
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: string_utl.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include "rc.h"
#include "string_utl.h"
#include "memory_manager.h"


/* 改行コードを取り除く */
/* MAC, DOS, UNIX 対応  */
/* 漢字コード未対応     */
void
clean_cr_string(char *s)
{
	int ii1 = 0;

	while((int)(s[ii1]) != '\0'){
		if( ((int)(s[ii1]) == 0x0a)
		  ||((int)(s[ii1]) == 0x0d) ){
			s[ii1] = (char)'\0';
		}
		ii1++;
	}
}


RC
init_string_array(STRING_ARRAY *stra)
{
	stra->size = 0;
	stra->alloc_size = 0;
	stra->str = NULL;

	return(NORMAL_RC);
}


RC
free_string_array(STRING_ARRAY *stra)
{
	if(stra->alloc_size > 0){
		if(stra->str == NULL) return(ARG_ERROR_RC);
		RC_TRY( mm_free(stra->str) );
	}
	stra->size = 0;
	stra->alloc_size = 0;
	stra->str = NULL;

	return(NORMAL_RC);
}


RC
add_string_array(STRING_ARRAY *stra, const char *string)
{
	int length;
	
	if(string == NULL) return(ARG_ERROR_RC);
	if(stra->size > stra->alloc_size) return(ARG_ERROR_RC);
	if(stra->size == stra->alloc_size){
		(stra->alloc_size) += 64;
		stra->str = mm_realloc(stra->str, (stra->alloc_size) * sizeof(char *));
		RC_NULL_CHK( stra->str );
	}
	length = strlen(string) + 1;
	if(length <= 0) return(UNKNOWN_ERROR_RC);
	RC_NULL_CHK( stra->str[stra->size] = mm_alloc(length * sizeof(char)) );
	strncpy(stra->str[stra->size], string, length);
	stra->str[stra->size][length-1] = (char)'\0';  /* 念の為 */
	(stra->size)++;

	return(NORMAL_RC);
}


RC
clean_string_array(STRING_ARRAY *stra)
{
	void **ptr;
	int ii1;

	if(stra->size > stra->alloc_size) return(ARG_ERROR_RC);
	if(stra->size < 0) return(ARG_ERROR_RC);

	if(stra->size > 0){
		RC_NULL_CHK( stra->str = mm_realloc(stra->str,
		                                    stra->size*sizeof(char *)) );
		stra->alloc_size = stra->size;
	}else{
		if(stra->alloc_size > 0) RC_TRY( mm_free(stra->str) );
		stra->str = NULL;
		stra->size = 0;
		stra->alloc_size = 0;
	}

	RC_NULL_CHK( ptr = mm_alloc_tmp((stra->size + 1)*sizeof(void *)) );
	for(ii1=0; ii1<stra->size; ii1++){
		ptr[ii1] = (void *)(stra->str[ii1]);
	}
	ptr[stra->size] = (void *)stra->str;
	RC_TRY( mm_compaction_arr(stra->size + 1, ptr) );
	stra->str = (char **)ptr[stra->size];
	for(ii1=0; ii1<stra->size; ii1++){
		stra->str[ii1] = (char *)ptr[ii1];
	}
	RC_TRY( mm_free(ptr) );

	return(NORMAL_RC);
}


/* str[size] の中で，文字列 s に一致するものを探し，そのインデックスを返却，*/
/* 見つからない場合は -1 を返却 */
int
search_string_array(int size, char **str, const char *s)
{
	int ii1;
	int ret = -1;

	for(ii1=0; ii1<size; ii1++){
		if(case_strcmp(str[ii1], s) == 0){
			ret = ii1;
			break;
		}
	}

	return(ret);
}


void
print_string_array(FILE *fp, STRING_ARRAY stra)
{
	int ii1;

	fprintf(fp, "size: %5d\n", stra.size);
	fprintf(fp, "alloc_size: %5d\n", stra.alloc_size);
	for(ii1=0;ii1<stra.size;ii1++){
		fprintf(fp, "string[%5d]: \"%s\"\n", ii1, stra.str[ii1]);
	}
}


RC
split_string_delim(STRING_ARRAY *stra, const char *string, const char *delim)
{
	int ii1, ii2, ii3;
	int str_size, delim_flag;
	char *s_tmp;

	if(string == NULL) return(ARG_ERROR_RC);
	if(delim == NULL) return(ARG_ERROR_RC);
	str_size = (strlen(string)+1);
	RC_NULL_CHK( s_tmp = mm_alloc(str_size*sizeof(char)) );
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
	RC_TRY( mm_free(s_tmp) );

	return(NORMAL_RC);
}


RC
strip_string_char(char *string, const char *strip)
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


RC
print_backspace(int num)
{
	if(num < 0) return(NORMAL_RC);
	for(; num; num--) if(fputc('\b', stderr) == EOF) return(WRITE_ERROR_RC);
	return(NORMAL_RC);
}


int
is_blank(const char *s)
{
	while(*s != (char)'\0'){
		if(*s != (char)' ') return(0);
		s++;
	}

	return(1);
}


int
case_strcmp(const char *s1, const char *s2)
{
	int ii1 = 0;
	int up1, up2;

	while(1){
		if( (s1[ii1] == (char)'\0')&&(s2[ii1] == (char)'\0') ) return(0);
		up1 = toupper((int)(s1[ii1]));
		up2 = toupper((int)(s2[ii1]));
		if(up1 != up2) return( up1 - up2 );
		ii1++;
	}

	return(0);
}


int
case_strncmp(const char *s1, const char *s2, size_t n)
{
	size_t ii1;
	int up1, up2;

	for(ii1=0; ii1<n; ii1++){
		if( (s1[ii1] == (char)'\0')&&(s2[ii1] == (char)'\0') ) return(0);
		up1 = toupper((int)(s1[ii1]));
		up2 = toupper((int)(s2[ii1]));
		if(up1 != up2) return( up1 - up2 );
	}

	return(0);
}


RC
str2double(const char *buf, double *val, int *p_flag)
{
	int sign = 1;
	int start_flag = 0;
	int stop_flag = 0;
	double point_val = 1.0;
	int point_flag = 0;
	int exp_flag = 0;
	int ii1 = 0;

	*val = 0.0;
	while(buf[ii1] != (char)'\0'){
		switch((int)(buf[ii1])){
		case ' ':
			if(start_flag) stop_flag = 1;
			break;
		case '+':
			if(start_flag){
				exp_flag = 1;
				stop_flag = 1;
			}
			break;
		case '-':
			if(start_flag){
				exp_flag = 1;
				stop_flag = 1;
			}else{
				sign *= -1;
			}
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			if(point_flag){
				point_val *= 0.1;
				*val += point_val * (double)((int)(buf[ii1]) - '0');
			}else{
				if(*val >= DBL_MAX/10.0001) return(OVERFLOW_ERROR_RC);
				*val *= 10.0;
				if(*val >= DBL_MAX - 2.0*((int)(buf[ii1]) - '0')){
					return(OVERFLOW_ERROR_RC);
				}
				*val += (double)((int)(buf[ii1]) - '0');
			}
			break;
		case '.':
			if(point_flag) return(CONVERT_ERROR_RC);
			point_flag = 1;
			break;
		case 'E':
		case 'e':
		case 'D':
		case 'd':
			ii1++;
			exp_flag = 1;
			stop_flag = 1;
			break;
		default:
			return(CONVERT_ERROR_RC);
		}
		if(stop_flag) break;
		ii1++;
	}
	if(start_flag == 0) return(CONVERT_ERROR_RC);
	if(p_flag != NULL) *p_flag = point_flag;

	if(exp_flag){
		int iexp;
		RC rc;

		rc = str2int(&(buf[ii1]), &iexp);
		if(rc != NORMAL_RC) return(rc);
		if(iexp > 0){
			for(ii1=0; ii1<iexp; ii1++){
				if(*val >= DBL_MAX/10.0001) return(OVERFLOW_ERROR_RC);
				*val *= 10.0;
			}
		}
		if(iexp < 0){
			for(ii1=0; ii1<(-iexp); ii1++){
				*val *= 0.1;
			}
		}
	}
	*val *= sign;

	return(NORMAL_RC);
}


RC
str2int(const char *buf, int *val)
{
	int sign = 1;
	int start_flag = 0;
	int stop_flag = 0;
	int ii1 = 0;

	*val = 0;
	while(buf[ii1] != (char)'\0'){
		switch((int)(buf[ii1])){
		case ' ':
			if(start_flag) stop_flag = 1;
			break;
		case '+':
			if(start_flag) return(CONVERT_ERROR_RC);
			break;
		case '-':
			if(start_flag){
				return(CONVERT_ERROR_RC);
			}else{
				sign *= -1;
			}
			break;
		case '0':
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
			start_flag = 1;
			if(*val > INT_MAX/10) return(OVERFLOW_ERROR_RC);
			*val *= 10;
			if(*val > INT_MAX - ((int)(buf[ii1]) - '0')){
				return(OVERFLOW_ERROR_RC);
			}
			*val += (int)(buf[ii1]) - '0';
			break;
		default:
			return(CONVERT_ERROR_RC);
		}
		if(stop_flag) break;
		ii1++;
	}
	if(start_flag == 0) return(CONVERT_ERROR_RC);

	*val *= sign;

	return(NORMAL_RC);
}


/* strncmp 文字列 buf, key の比較を行う．                         */
/* ただし、比較は最大 n 文字まで．                                */
/* buf の先頭の空白は無視し、無視した空白の数を skip_count に代入 */
int
space_strncmp(const char *buf, const char *key, size_t n, int *skip_count)
{
	int pos = 0;

	while((int)(buf[pos]) == ' '){
		pos++;
	}
	if(skip_count != NULL) *skip_count = pos;

/*	return(strncmp(&(buf[pos]), key, n));*/
	return(case_strncmp(&(buf[pos]), key, n));
}


/* ファイル名 fname[] の拡張子を除いて base[base_size] にコピー */
RC
set_base_name(const char fname[], char base[], int base_size)
{
	int ii1;
	int ptr, len;

	len = strlen(fname);
	ptr = len - 1;
	while(ptr >= 0){
		if(fname[ptr] == (char)'.') break;
		ptr--;
	}
	if(ptr < 0) ptr = len;
	if(ptr >= base_size) ptr = base_size - 1;

	for(ii1=0; ii1<ptr; ii1++){
		base[ii1] = fname[ii1];
	}
	base[ptr] = (char)'\0';

	return(NORMAL_RC);
}


RC
copy_file(const char fname_src[], const char fname_dest[])
{
	FILE *fp_src;
	FILE *fp_dest;
	char buf[512];

	RC_NULL_CHK( fp_src = fopen(fname_src, "rb") );
	RC_NULL_CHK( fp_dest = fopen(fname_dest, "wb") );

	while(1){
		size_t read_size = fread(buf, 1, sizeof(buf), fp_src);

		if(read_size == 0) break;
		RC_TRY( rc_fwrite(buf, 1, read_size, fp_dest) );
		if(read_size != sizeof(buf)) break;
	}

	fclose(fp_src);
	fclose(fp_dest);

	return(NORMAL_RC);
}


