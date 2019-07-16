/*********************************************************************
 * string_utl.h
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Ryu SASAOKA> <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: string_utl.h,v 1.7 2004/02/05 07:15:02 sasaoka Exp $ */

#ifndef STRING_UTL_H
#define STRING_UTL_H

#include "rc.h"

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
RC logprintf(char *f_name, char *format, ...);
RC init_string_array(STRING_ARRAY *stra);
RC free_string_array(STRING_ARRAY *stra);
RC add_string_array(STRING_ARRAY *stra, const char *string);
void print_string_array(FILE *fp, STRING_ARRAY stra);
RC split_string_delim(STRING_ARRAY *stra, const char *string,const char *delim);
RC strip_string_char(char *string, const char *strip);
RC print_backspace(int num);
RC print_progress(const char *format, ...);
RC print_progress_bar(double percent, const char *format, ...);

/* idc.c */
RC fgetidc(FILE *fp, const char coord[], IDC data[]);
RC sgetidc(const char *s, const char coord[], IDC data[]);

#endif /* STRING_UTL_H */


