/*********************************************************************
 * nst_blk_llio.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nst_blk_llio.c,v 1.13 2003/07/23 04:06:31 sasaoka Exp $ */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include "rc.h"
#include "string_utl.h"
#include "nst_component.h"

#define BUFFER_SIZE  (512)
#define NST_TOK_SIZE (24)

/* local functions */
static int inc_bufs_point(int point);
static int dec_bufs_point(int point);
static void clean_cr_string(char *s);
static FILE *include_fopen(const char *buf, int skip);
static int space_strncmp(const char *buf, const char *key,
                         size_t n, int *skip_count);
static int nsttok(const char *buf, int offset, char *tok,
                  int max_size, int *eol_flag);
static void supple_E(const char *src, char *dest);
static void clean_n_string(char *str, int n);
static char *d2s8(double val, char *buf);


RC init_nst_file(FILE *fp, NST_FILE *nst_file)
{
	int ii1;

	rewind(fp);
	nst_file->parent = fp;
	nst_file->child = NULL;

	for(ii1=0;ii1<NST_FILE_BUFS;ii1++){
		nst_file->bufs[ii1] = (char *)malloc(BUFFER_SIZE * sizeof(char));
		if(nst_file->bufs[ii1] == NULL){
			fprintf(stderr, "malloc < ");
			return(ALLOC_ERROR_RC);
		}
		nst_file->bufs[ii1][0] = (char)'\0';
	}
	nst_file->r_point = nst_file->w_point = 0;

	return(NORMAL_RC);
}


RC term_nst_file(NST_FILE *nst_file)
{
	int ii1;

	if(nst_file->child != NULL){
		fclose(nst_file->child);
		nst_file->child = NULL;
	}

	for(ii1=0;ii1<NST_FILE_BUFS;ii1++){
		free(nst_file->bufs[ii1]);
	}

	return(NORMAL_RC);
}


static int inc_bufs_point(int point)
{
	point++;
	if(point >= NST_FILE_BUFS){
		point = 0;
	}

	return(point);
}


static int dec_bufs_point(int point)
{
	point--;
	if(point < 0){
		point = NST_FILE_BUFS - 1;
	}

	return(point);
}


RC nst_allocate_tok(char **toks, int size)
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);

	toks[0] = malloc(size * NST_TOK_SIZE * sizeof(char));
	if(toks[0] == NULL){
		fprintf(stderr, "malloc < ");
		return(ALLOC_ERROR_RC);
	}

	for(ii1=1; ii1<size; ii1++){
		toks[ii1] = toks[0] + NST_TOK_SIZE * ii1;
	}

	return(NORMAL_RC);
}


RC nst_free_tok(char **toks, int size)
{
	free((void *)toks[0]);

	return(NORMAL_RC);
}


void nst_print_tok(char **toks, int size)
{
	int ii1;

	for(ii1=0; ii1<size; ii1++){
		fprintf(stdout, "[%d][%s]\n", ii1, toks[ii1]);
	}
}


RC nst_seek_begin_bulk(NST_FILE *nst_file)
{
	char *buf;

	while(1){
		RC_TRY( nst_gets(&buf, nst_file) );
		if(nst_chk_begin_bulk(buf)) break;
	}

	return(NORMAL_RC);
}


RC nst_seek_cend(NST_FILE *nst_file)
{
	char *buf;

	while(1){
		RC_TRY( nst_gets(&buf, nst_file) );
		if(nst_chk_cend(buf)) break;
	}

	return(NORMAL_RC);
}


/* NASTRAN bulk data file $B$N(B INCLUDE $BJ8$r9MN8$7$F#19TFI$_9~$`(B */
/* s $B$KFI$_9~$s$@%P%C%U%!$N%]%$%s%?$rBeF~$9$k(B */
RC nst_gets(char **s, NST_FILE *nst_file)
{
	int skip_count;

	if( (nst_file->r_point >= NST_FILE_BUFS)
	  ||(nst_file->w_point >= NST_FILE_BUFS)
	  ||(nst_file->r_point < 0)||(nst_file->w_point < 0) ){
		return(ARG_ERROR_RC);
	}

	if(nst_file->r_point != nst_file->w_point){
		*s = nst_file->bufs[nst_file->r_point];
		nst_file->r_point = inc_bufs_point(nst_file->r_point);
		return(NORMAL_RC);
	}

	*s = nst_file->bufs[nst_file->w_point];

	/* child $B$+$i$NFI$_9~$_(B */
	if(nst_file->child != NULL){	
		if(fgets(*s, BUFFER_SIZE, nst_file->child) == NULL){
			/* child $B$+$i$NFI$_9~$_$K<:GT(B */
			if(feof(nst_file->child)){
				/* $B%U%!%$%k=*C<$KC#$7$?!#(B */
				fclose(nst_file->child);
				nst_file->child = NULL;
				return(nst_gets(s, nst_file));
			}
			/* $BFI$_9~$_%(%i!<(B */
			return(READ_ERROR_RC);
		}
		(*s)[BUFFER_SIZE - 1] = '\0';
		clean_cr_string(*s);
		nst_file->r_point = inc_bufs_point(nst_file->r_point);
		nst_file->w_point = inc_bufs_point(nst_file->w_point);
		return(NORMAL_RC);
	}

	/* parent $B$+$i$NFI$_9~$_(B */
	if(fgets(*s, BUFFER_SIZE, nst_file->parent) == NULL){
		if(feof(nst_file->parent)){
			return(END_RC);
		}
		return(READ_ERROR_RC);
	}
	(*s)[BUFFER_SIZE - 1] = '\0';
	clean_cr_string(*s);

	if(space_strncmp(*s, "INCLUDE", 7, &skip_count) == 0){
		/* INCLUDE $B9T$N=hM}(B */
		nst_file->child = include_fopen(*s, 7+skip_count); 
		if(nst_file->child == NULL){
			fprintf(stderr, "include_fopen < ");
			return(OPEN_ERROR_RC);
		}
		return(nst_gets(s, nst_file));
	}

	nst_file->r_point = inc_bufs_point(nst_file->r_point);
	nst_file->w_point = inc_bufs_point(nst_file->w_point);
	return(NORMAL_RC);
}


/* $B%U%!%$%k%]%$%s%?$r!J2>A[E*$K!K0l9TJ,La$9(B */
RC nst_ungets(NST_FILE *nst_file)
{
	nst_file->r_point = dec_bufs_point(nst_file->r_point);

	return(NORMAL_RC);
}


/* $B2~9T%3!<%I$r<h$j=|$/(B */
/* MAC, DOS, UNIX $BBP1~(B */
/* $B4A;z%3!<%IL$BP1~(B */
static void clean_cr_string(char *s)
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


static FILE *include_fopen(const char *buf, int skip)
{
	int index = skip;
	int wpos = 0;
	int q_flag = 0;
	char name[BUFFER_SIZE];

	while((int)(buf[index]) != '\0'){
		if((int)(buf[index]) == '\''){
			if(q_flag){
				break;
			}else{
				q_flag = 1;
			}
		}else{
			if(((int)(buf[index]) != ' ')||(q_flag)){
				name[wpos] = buf[index];
				wpos++;
			}
		}
		index ++;
	}
	name[wpos] = (char)'\0';

	return(fopen(name, "r"));
}


/* strncmp $BJ8;zNs(B buf, key $B$NHf3S$r9T$&(B */
/* $B$?$@$7!"Hf3S$O:GBg(B n $BJ8;z$^$G!#(B*/
/* buf $B$N@hF,$N6uGr$OL5;k$7!"L5;k$7$?6uGr$N?t$r(B skip_count $B$KBeF~(B */
static int space_strncmp(const char *buf, const char *key,
                         size_t n, int *skip_count)
{
	int pos = 0;

	while((int)(buf[pos]) == ' '){
		pos++;
	}
	if(skip_count != NULL) *skip_count = pos;

	return(strncmp(&(buf[pos]), key, n));
}



/* buf $B$,(B BEGIN BULK $B9T$J$i(B 1$B!"$=$&$G$J$$$J$i(B 0 $B$rJV5Q(B */
int nst_chk_begin_bulk(const char *buf)
{
	int skip_count0;

	if(space_strncmp(buf, "BEGIN", 5, &skip_count0) == 0){
		if(space_strncmp(buf + skip_count0 + 5,
		                 "BULK", 4, NULL) == 0){
			return(1);
		}
	}

	return(0);
}


/* buf $B$,(B BEGIN BULK $B9T$J$i(B 1$B!"$=$&$G$J$$$J$i(B 0 $B$rJV5Q(B */
int nst_chk_cend(const char *buf)
{
	if(space_strncmp(buf, "CEND", 4, NULL) == 0){
		return(1);
	}

	return(0);
}


RC nst_blk_get_tok(NST_FILE *nst_file,
                   int num_keys, NST_BULK_KEY *keys, int *hit_key,
                   int max_toks, char **toks, int *num_toks)
{
	int ii1;
	RC rc;
	char *buf;
	int offset;
	int skip_count;
	int long_flag = 0;
	int max_line_toks = 0;
	int max_tok_size = 0;
	int eol_flag = 0;

	if((num_keys <= 0)||(max_toks <= 0)){
		return(ARG_ERROR_RC);
	}

	while(1){  /* keys $B$H0lCW$9$k9T$,8=$l$k$^$G(B */
		rc = nst_gets(&buf, nst_file);
		if(rc != NORMAL_RC){
			return(rc);
		}
		offset = nsttok(buf, 0, toks[0], 8, &eol_flag);
		if(toks[0][0] == (char)'$') continue;   /* $B%3%a%s%H9T$rL5;k(B */
		if(eol_flag) continue;  /* $B%G!<%?$N$J$$9T$rL5;k(B */

		*hit_key = -1;
		for(ii1=0;ii1<num_keys;ii1++){
			if(space_strncmp(toks[0], keys[ii1].key,
			                 keys[ii1].key_size, &skip_count) == 0){
			     if((int)(toks[0][skip_count + keys[ii1].key_size]) == ' '){
					*hit_key = ii1;
					long_flag = 0;
					break;
				}else if((int)(toks[0][skip_count
				                       + keys[ii1].key_size]) == '*'){
					*hit_key = ii1;
					long_flag = 1;
					break;
				}
			}
		}
		if(*hit_key >= 0) break;
	}

	*num_toks = 1;

	if(long_flag){
		max_line_toks = 4;
		max_tok_size = 16;
	}else{
		max_line_toks = 8;
		max_tok_size = 8;
	}
	while(1){
		for(ii1=0;ii1<max_line_toks;ii1++){
			if(*num_toks >= max_toks){
				return(OVERFLOW_ERROR_RC);
			}
			offset = nsttok(buf, offset, toks[*num_toks],
			                max_tok_size, &eol_flag);
			(*num_toks)++;
			/* if(eol_flag) break; */
		}

		/* $B5l%U%)!<%^%C%HMQ(B
		if(eol_flag) break;
		offset = nsttok(buf, offset, toks[*num_toks],
		                max_tok_size, &eol_flag);
		if( ((int)(toks[*num_toks][0]) != '+')
		  &&((int)(toks[*num_toks][0]) != '*')){
			break;
		}
		*/

		while(1){
			rc = nst_gets(&buf, nst_file);
			if(rc != NORMAL_RC){
				fprintf(stderr, "nst_gets < ");
				if(rc == END_RC) return(READ_ERROR_RC);
				return(rc);
			}
			if(*num_toks >= max_toks){
				return(OVERFLOW_ERROR_RC);
			}
			offset = nsttok(buf, 0, toks[*num_toks], 8, &eol_flag);
			if(toks[*num_toks][0] != (char)'$') break; /* $B%3%a%s%H9T$rL5;k(B */
		}
		if( (((int)(toks[*num_toks][0]) != '+')
		     &&((int)(toks[*num_toks][0]) != '*')
		     &&((int)(toks[*num_toks][0]) != ' ') ) || eol_flag){
			/* return(READ_ERROR_RC); $B5l%U%)!<%^%C%HMQ(B */
			nst_ungets(nst_file);
			break;
		}
		if((int)(toks[*num_toks][0]) == '*'){
			max_line_toks = 4;
			max_tok_size = 16;
		}else{
			max_line_toks = 8;
			max_tok_size = 8;
		}
	}

	return(NORMAL_RC);
}


/* $BCm0U(B: tok $B$O(B max_size+2 $BJ8;z0J>e3JG<$G$-$J$1$l$P$J$i$J$$(B */
/* $B!J(Btok $B$N:G8e$O0lJ8;z0J>e$N6uGr(B + '\0' $B$G$"$k$3$H$rJ]>Z$9$k$?$a!K(B */
static int nsttok(const char *buf, int offset, char *tok,
                  int max_size, int *eol_flag)
{
	int rp = offset;
	int wp = 0;
	int w_flag = 0;
	int r_count = 0;

	while((int)(buf[rp]) != '\0'){
		if((int)(buf[rp]) == ','){
			/* Free Field Format $BBP:v!J<c43LdBj$"$j!K(B */
			rp++;
			break;
		}
		if((w_flag != 0)||((int)buf[rp] != ' ')){
			w_flag = 1;
			tok[wp] = buf[rp];
			wp++;
		}
		rp++;
		r_count++;
		if(r_count >= max_size){
			break;
		}
	}

	for(;wp<=max_size;wp++){
		tok[wp] = (char)' ';
	}
	tok[max_size+1] = (char)'\0';

	if((int)buf[rp] == '\0'){
		*eol_flag = 1;
	}else{
		*eol_flag = 0;
	}

	return(rp);
}


double nst_tok2d(int tok_size, int tok_num, char **toks)
{
	char buf[BUFFER_SIZE];
	double ret;

	if(tok_num >= tok_size){
		return(FNSDEFAULT);
	}

	supple_E(toks[tok_num], buf);
	if(sscanf(buf, "%le", &ret) != 1){
		return(FNSDEFAULT);
	}

	return(ret);
}


RC nst_get_param(const char name[], const char format[],
                 NST_PARAM_ARRAY param, void *ptr1, void *ptr2)
{
	int ii1;
	int skip_count;
	int name_size = strlen(name);
	int index;
	char buf[BUFFER_SIZE];

	for(ii1=0; ii1<param.size; ii1++){
		if(space_strncmp(param.array[ii1].name, name,
		                 name_size, &skip_count) == 0){
			if( (param.array[ii1].name[name_size+skip_count] == ' ')
			  ||(param.array[ii1].name[name_size+skip_count] == '\0') ){
				break;
			}
		}
	}
	index = ii1;

	if(ptr1 != NULL){
		if(format[0] == 'i'){
			if(index >= param.size){
				*((int *)ptr1) = NSDEFAULT;
			}else{
				if(sscanf(param.array[index].val_1.s, "%d", (int *)ptr1) != 1){
					*((int *)ptr1) = NSDEFAULT;
				}
			}
		}else if(format[0] == 'd'){
			if(index >= param.size){
				*((double *)ptr1) = FNSDEFAULT;
			}else{
				supple_E(param.array[index].val_1.s, buf);
				if(sscanf(buf, "%le", (double *)ptr1) != 1){
					*((double *)ptr1) = FNSDEFAULT;
				}
			}
		}
	}

	if( (ptr2 != NULL)&&(format[0] != '\0') ){
		if(format[1] == 'i'){
			if(index >= param.size){
				*((int *)ptr2) = NSDEFAULT;
			}else{
				if(sscanf(param.array[index].val_2.s, "%d", (int *)ptr2) != 1){
					*((int *)ptr2) = NSDEFAULT;
				}
			}
		}else if(format[1] == 'd'){
			if(index >= param.size){
				*((double *)ptr2) = FNSDEFAULT;
			}else{
				supple_E(param.array[index].val_2.s, buf);
				if(sscanf(buf, "%le", (double *)ptr2) != 1){
					*((double *)ptr2) = FNSDEFAULT;
				}
			}
		}
	}

	return(NORMAL_RC);
}


/* supplement E of exponential expression, if it is missing */
static void supple_E(const char *src, char *dest)
{
	int rp = 0;
	int wp = 0;
	int num_flag = 0;
	int E_flag = 0;

	while((int)(src[rp]) != '\0'){
		if(('0' <= (int)(src[rp]))&&((int)(src[rp]) <= '9')){
			num_flag = 1;
		}
		if(((int)(src[rp]) == 'E')||((int)(src[rp]) == 'e')
		 ||((int)(src[rp]) == 'D')||((int)(src[rp]) == 'd')){
			if((num_flag)&&(!E_flag)){
				E_flag = 1;
				dest[wp] = (char)'E';
				rp++;
				wp++;
				continue;
			}
		}
		if(((int)(src[rp]) == '+')||((int)(src[rp]) == '-')){
			if((num_flag)&&(!E_flag)){
				dest[wp] = (char)'E';
				wp++;
				E_flag = 1;
			}
		}
		dest[wp] = src[rp];
		rp++;
		wp++;
	}
	dest[wp] = (char)'\0';
}


int nst_tok2i(int tok_size, int tok_num, char **toks)
{
	int ret;

	if(tok_num >= tok_size){
		return(NSDEFAULT);
	}

	if(sscanf(toks[tok_num], "%d", &ret) != 1){
		return(NSDEFAULT);
	}

	return(ret);
}


int nst_tok_strncmp(int tok_size, int tok_num, char **toks,
                                 const char *str, int size)
{
	int ret;

	if(tok_num >= tok_size){
		return(1);
	}
	ret = strncmp(toks[tok_num], str, size);
	if(ret) return(ret);

	if(toks[tok_num][size] == (char)' ') return(0);

	return(-1);
}


RC nst_print(FILE *stream, const char *format, const IDC *data)
{
	int ii1 = 0;
	int ii2;
	int index = 0;
	char buf[BUFFER_SIZE];

	while(format[ii1]){
		switch((int)format[ii1]){
		case 'E':
			/* $BJ8;zNsA4BN(B(Entire string) */
			fprintf(stream, "%s", data[index].c);
			index++;
			break;

		case 'i':
			/* 8 $B%+%i%`@0?t(B */
			if(data[index].i == (int)NSDEFAULT){
				fprintf(stream,"        ");
			}else{
				fprintf(stream, "%8d", data[index].i);
			}
			index++;
			break;

		case 'I':
			/* 16 $B%+%i%`@0?t(B */
			if(data[index].i == (int)NSDEFAULT){
				fprintf(stream,"                ");
			}else{
				fprintf(stream, "%16d", data[index].i);
			}
			index++;
			break;

		case 'd':
			/* 8 $B%+%i%`<B?t(B */
			if(data[index].d <= CR_FNSDEFAULT){
				fprintf(stream, "        ");
			}else{
				fprintf(stream,"%s", d2s8(data[index].d, buf));
			}
			index++;
			break;

		case 'D':
			/* 16 $B%+%i%`<B?t(B */
			if(data[index].d <= CR_FNSDEFAULT){
				fprintf(stream,"                ");
			}else if(fabs(data[index].d) < 2.0E-99){
				fprintf(stream, " 0.0000000000000");
			}else{
				fprintf(stream, "%16.8E", data[index].d);
			}
			index++;
			break;

		case 'G':
			/* 16 $B%+%i%`<B?t(B($BG\@:EYI=<((B) */
			if(data[index].d <= CR_FNSDEFAULT){
				fprintf(stream,"                ");
			}else if(fabs(data[index].d) < 2.0E-99){
				sprintf(buf, "%16.8E", 0.0);
			}else{
				sprintf(buf, "%16.8E", data[index].d);
			}
			for(ii2=0; ii2<16; ii2++){
				if(buf[ii2] == (char)'E') buf[ii2] = (char)'D';
			}
			fprintf(stream, "%s", buf);
			index++;
			break;

		case 'c':
			/* 8 $B%+%i%`J8;zNs(B */
			strcpy(buf, data[index].c);
			index++;
			clean_n_string(buf, 8);
			fprintf(stream, "%s", buf);
			break;

		case 'C':
			/* 16 $B%+%i%`J8;zNs(B */
			strcpy(buf, data[index].c);
			index++;
			clean_n_string(buf, 16);
			fprintf(stream, "%s", buf);
			break;

		case 's':
			/* 8 $B%+%i%`J8;zNs(B */
			strcpy(buf, data[index].s);
			index++;
			clean_n_string(buf, 8);
			fprintf(stream, "%s", buf);
			break;

		case 'S':
			/* 16 $B%+%i%`J8;zNs(B */
			strcpy(buf, data[index].s);
			index++;
			clean_n_string(buf, 16);
			fprintf(stream, "%s", buf);
			break;

		case '-':
			/* 8 $B%+%i%`6uGr(B */
			fprintf(stream, "        ");
			break;

		case '=':
			/* 16 $B%+%i%`6uGr(B */
			fprintf(stream,"                ");
			break;

		case 'r':
			fprintf(stream, "\n");
			break;

		default:
			return(ARG_ERROR_RC);
		}

		ii1++;
	}

	return(NORMAL_RC);
}


/* convert to just n character string  */
static void clean_n_string(char *str, int n)
{
	int rp,wp;

	/* delete '\n' character */
	str[n] = '\0';
	clean_cr_string(str);

	/* delete ' ' character */
	wp = 0;
	rp = 0;
	while(str[rp] != '\0'){
		if(str[rp] != ' '){
			str[wp] = str[rp];
			wp++;
		}
		rp++;
	}
	for(;wp<n;wp++){
		str[wp] = ' ';
	}
}


/* double -> string(8char) */
/* double $B7?JQ?t(B val $B$r(B 8 $BJ8;z$NJ8;zNs$KJQ49$7!"(Bbuf $B$K=q$-9~$`(B  */
/* buf $B$NCM$,$=$N$^$^JV5Q$5$l$k(B */
static char *d2s8(double val, char *buf)
{
	double abs_val;
	int precision;

	buf[0] = '\0';
	abs_val = fabs(val);

	if(abs_val < (2.0E-99)){
		/* nearly zero */
		sprintf(buf,"0.000000");
	}else if(val >= (9.98E+99)){
		sprintf(buf,"9.99E+99");
	}else if(val <= (-9.8E+99)){
		sprintf(buf,"-9.9E+99");
	}else if(abs_val < 0.001){
		precision = 2;
		if(val > 0.0)precision++;
		sprintf(buf,"%9.*E",precision,val);
		if(strlen(buf) > 9){    /* For Windows VC++ */
			buf[7] = buf[8];
			buf[8] = buf[9];
			buf[9] = '\0';
		}
		if(buf[7] == '0'){
			buf[7] = buf[8];
			buf[8] = '\0';
		}else{
			precision = 1;
			if(val > 0.0)precision++;
			sprintf(buf,"%8.*E",precision,val);
			if(strlen(buf) > 8){    /* For Windows VC++ */
				buf[6] = buf[7];
				buf[7] = buf[8];
				buf[8] = '\0';
			}
		}
	}else if(abs_val < 999999.){
		if(abs_val < 9.99999){
			precision = 5;
		}else if(abs_val < 99.9999){
			precision = 4;
		}else if(abs_val < 999.999){
			precision = 3;
		}else if(abs_val < 9999.99){
			precision = 2;
		}else if(abs_val < 99999.9){
			precision = 1;
		}else{
			precision = 0;
		}
		if(val > 0.0)precision++;
		if(precision == 0){
			sprintf(buf,"%7.0f.",val);
		}else{
			sprintf(buf,"%8.*f",precision,val);
		}
	}else{
		precision = 2;
		if(val > 0.0)precision++;
		sprintf(buf,"%9.*E",precision,val);
		if(strlen(buf) > 9){    /* For Windows VC++ */
			buf[7] = buf[8];
			buf[8] = buf[9];
			buf[9] = '\0';
		}
		if(buf[7] == '0'){
			buf[7] = buf[8];
			buf[8] = '\0';
		}else{
			precision = 1;
			if(val > 0.0)precision++;
			sprintf(buf,"%8.*E",precision,val);
			if(strlen(buf) > 8){    /* For Windows VC++ */
				buf[6] = buf[7];
				buf[7] = buf[8];
				buf[8] = '\0';
			}
		}
	}
	
	buf[8] = '\0';

	return(buf);
}











#if 0


int main(void)
{
	FILE *fp;
	NST_FILE nst_file;
	int toks_size = 64;
	char *toks[64];
	int num_keys = 2;
	NST_BULK_KEY keys[] = {{"SPC1", 4}, {"SPCADD", 6}};

	fp = fopen("hoge.dat", "r");
	if(fp == NULL){
		fprintf(stderr, "Error\n");
		return(0);
	}
	if(nst_allocate_tok(toks, toks_size) != NORMAL_RC){
		fprintf(stderr, "nst_allocate_tok : \n");
		return(0);
	}
	init_nst_file(fp, &nst_file);
	if(nst_seek_begin_bulk(&nst_file) != NORMAL_RC){
		fprintf(stderr, "seek_begin_bulk\n");
		return(0);
	}
	while(1){
		RC rc;
		int hit_key;
		int num_toks;
		int ii1;

		rc = nst_get_blk(&nst_file, num_keys, keys, &hit_key,
		                 toks_size, toks, &num_toks);
		if(rc != NORMAL_RC){
			if(rc == END_RC) break;
			fprintf(stderr, "nst_get_blk : \n");
			return(0);
		}

		for(ii1=0;ii1<num_toks;ii1++){
			fprintf(stdout, "%2d/%2d:[%s]\n", ii1, num_toks, toks[ii1]);
		}
	}
	nst_free_tok(toks, toks_size);

	return(0);
}


#endif /* 0 */

/*
memo:
*/



