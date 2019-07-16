/*********************************************************************
 * string_utl.h
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: string_utl.h 339 2004-07-14 05:17:44Z sasaoka $ */

#ifndef STRING_UTL_H
#define STRING_UTL_H

#include <stdlib.h>
#include "rc.h"
#include "memory_manager.h"


/* string_utl.c */
typedef struct{
	int size;
	int alloc_size;
	char **str;
} STRING_ARRAY;


/* idc.c */
#define MAX_IDC_S_SIZE 64

typedef struct{
	int size;
	long double *v;
} LONG_DOUBLE_ARRAY;

typedef union{
	int i;
	long int li;
	double d;
	long double ld;
	LONG_DOUBLE_ARRAY ldarr;
	char *c;
	char s[MAX_IDC_S_SIZE+1];
	char a;
}IDC;


/* string_utl.c */
void clean_cr_string(char *s);
RC init_string_array(STRING_ARRAY *stra);
RC free_string_array(STRING_ARRAY *stra);
RC add_string_array(STRING_ARRAY *stra, const char *string);
RC clean_string_array(STRING_ARRAY *stra);
void print_string_array(FILE *fp, STRING_ARRAY stra);

RC split_string_delim(STRING_ARRAY *stra, const char *string,const char *delim);
RC strip_string_char(char *string, const char *strip);
RC print_backspace(int num);
int is_blank(const char *s);
int case_strcmp(const char *s1, const char *s2);
int case_strncmp(const char *s1, const char *s2, size_t n);
RC str2double(const char *buf, double *val, int *p_flag);
RC str2int(const char *buf, int *val);
int space_strncmp(const char *buf, const char *key, size_t n, int *skip_count);
RC set_base_name(const char fname[], char base[], int base_size);


/* idc.c */
RC fgetidc(FILE *fp, const char coord[], IDC data[]);
RC sgetidc(const char *s, const char coord[], IDC data[]);

#endif /* STRING_UTL_H */


