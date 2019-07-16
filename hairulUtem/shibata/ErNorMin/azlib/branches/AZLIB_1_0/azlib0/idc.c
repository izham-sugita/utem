/*********************************************************************
 * idc.c
 *
 * Copyright (C) 2002 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Mitsunori HIROSE>
 *
 * Thanks to following contributors.
 *  <Ryu SASAOKA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: idc.c,v 1.6 2003/07/22 10:21:27 adachi Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "rc.h"
#include "string_utl.h"

#define BUFFER_SIZE 1024

static int idc_get_c(const char buf[], char tmps[], int ptr);
static int idc_get_a(const char buf[], char tmps[], int ptr);
static int idc_get_d(const char buf[], char tmps[], int ptr);
static int idc_get_i(const char buf[], char tmps[], int ptr);

#ifdef DEBUG_MAIN
int main(void);
int main(void)
{
	IDC data[10];
	char *str = " XABCDEFG 1.0D-10-11-12.0 1,2  3.0XYZUVW 2.54";

	printf("%s\n",str);
	sgetidc(str,"acdddiiisi",data);

	printf("\n");
	printf("data[0].a = %c\n",data[0].a);
	printf("data[1].c = %s\n",data[1].c);
	printf("data[2].d = %e\n",data[2].d);
	printf("data[3].d = %e\n",data[3].d);
	printf("data[4].d = %e\n",data[4].d);
	printf("data[5].i = %d\n",data[5].i);
	printf("data[6].i = %d\n",data[6].i);
	printf("data[7].i = %d\n",data[7].i);
	printf("data[8].s = %s\n",data[8].s);
	printf("data[9].i = %d\n",data[9].i);

	return(0);
}
#endif /* DEBUG_MAIN */


RC fgetidc(FILE *fp,const char coord[], IDC data[])
{
	char buf[BUFFER_SIZE];
	
	if(fgets(buf,sizeof(buf),fp) == NULL){
		if(feof(fp)){
			return(END_RC);
		}
		return(READ_ERROR_RC);
	}

	return(sgetidc(buf,coord,data));
}

RC sgetidc(const char buf[], const char coord[], IDC data[])
{
	char tmps[BUFFER_SIZE];
	int ii = 0;
	int ptr = 0;
	double d;
	int i;
	char *c;
	int len;

	while(coord[ii] != '\0'){
		if(coord[ii] == 'd'){
			ptr = idc_get_d(buf,tmps,ptr);
#ifdef DEBUG_MAIN
printf(" [d|%s]",tmps);
#endif
			if(sscanf(tmps,"%le",&d) != 1){
				return(CONVERT_ERROR_RC);
			}
			data[ii].d = d;
		}else if(coord[ii] == 'i'){
			ptr = idc_get_i(buf,tmps,ptr);
#ifdef DEBUG_MAIN
printf(" [i|%s]",tmps);
#endif
			if(sscanf(tmps,"%d",&i) != 1){
				return(CONVERT_ERROR_RC);
			}
			data[ii].i = i;
		}else if((coord[ii] == 'c')||(coord[ii] == 's')){
			ptr = idc_get_c(buf,tmps,ptr);
			len = strlen(tmps);
			if(len <= 0){
				return(CONVERT_ERROR_RC);
			}
			if(coord[ii] == 'c'){
#ifdef DEBUG_MAIN
printf(" [c|%s]",tmps);
#endif
				if((c = malloc(len+1)) == NULL){
					return(ALLOC_ERROR_RC);
				}
				strcpy(c,tmps);
				data[ii].c = c;
			}else{   /* (coord[ii] == 's') */
				tmps[MAX_IDC_S_SIZE] = '\0';
#ifdef DEBUG_MAIN
printf(" [s|%s]",tmps);
#endif
				strcpy(data[ii].s, tmps);
			}
		}else if(coord[ii] == 'a'){
			ptr = idc_get_a(buf,tmps,ptr);
#ifdef DEBUG_MAIN
printf(" [a|%s]",tmps);
#endif
			if(tmps[0] == '\0'){
				return(CONVERT_ERROR_RC);
			}
			data[ii].a = tmps[0];
		}else{
			return(ARG_ERROR_RC);
		}
		ii++;
	}
	return(NORMAL_RC);
}


/* buf[ptr] 以降で空白文字等で区切られた文字列をtmpsに切り出す */
static int idc_get_c(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* 文字列が現れたら1 */
	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;  /* 文字列(buf)の終端に達した */
		}
		if((buf[ptr] == ' ')||(buf[ptr] == '\t')||(buf[ptr] == ',')){
			/* 区切りになる文字が現れた */
			if(vflag == 1){
				/* すでに文字列を読んだのでループを出る */
				break;
			}else{
				/* まだ文字列を読んでないのでとりあえず先に進む */
				ptr++;
			}
		}else{
			/* 一般的な文字が現れた */
			vflag = 1;
			tmps[ii] = buf[ptr]; /* １文字コピー */
			ptr++;
			ii++;
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}
	tmps[ii] = '\0';

	return(ptr);
}

/* buf[ptr] 以降で空白文字等以外の１文字をtmps[0]に切り出す */
static int idc_get_a(const char buf[], char tmps[], int ptr)
{
	int ii = 0;

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;  /* 文字列(buf)の終端に達した */
		}
		if((buf[ptr] == ' ')||(buf[ptr] == '\t')||(buf[ptr] == ',')){
			/* 空白文字が現れた */
			/* 先に進む */
			ptr++;
		}else{
			/* 一般的な文字が現れた */
			tmps[ii] = buf[ptr]; /* １文字コピー */
			ptr++;
			ii++;
			break;
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}
	tmps[ii] = '\0';

	return(ptr);
}



/* buf[ptr] 以降で数値(らしきもの)をtmpsに切り出す */
static int idc_get_d(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* 数値らしきものが現れたら1 */
	int sflag = 0;  /* 符号が現れたら1 */
	int eflag = 0;  /* 指数部らしきものが現れたら1 */
	int pflag = 0;  /* 小数点らしきものが現れたら1 */

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;   /* 文字列の終端に達した */
		}
		if((buf[ptr] == ',')||(buf[ptr] == ' ')||(buf[ptr] == '\t')){
			/* 区切りになる文字が現れた */
			if(vflag == 1){
				/* すでに数値らしきものを読んだのでループを出る */
				break;
			}else{
				/* まだ数値らしきものを読んでないのでとりあえず先に進む */
				ptr++;
			}
		}else if(('0' <= buf[ptr])&&(buf[ptr] <= '9')){
			/* 数字が現れた */
			vflag = 1;
			sflag = 1; /* 符号も現れたことにする */
			tmps[ii] = buf[ptr]; /* １文字コピー */
			ptr++;
			ii++;
		}else if((buf[ptr] == '-')||(buf[ptr] == '+')){
			/* 符号が現れた */
			if(sflag == 1){
				/* 既に符号が現れている */
				break;
			}else{
				/* 符号が現れたのははじめて */
				sflag = 1;
				tmps[ii] = buf[ptr]; /* １文字コピー */
				ptr++;
				ii++;
			}
		}else if((buf[ptr] == 'D')||(buf[ptr] == 'E')
		       ||(buf[ptr] == 'd')||(buf[ptr] == 'e')){
			/* D とか E が現れた */
			if(vflag == 0){
				/* まだ数字を読んでいないので無視して次に進む */
				ptr++;
			}else{
				if(eflag == 0){
					eflag = 1; /* これから指数部を読む */
					sflag = 0; /* 符号のフラグをりセット */
					pflag = 1; /* 小数点のフラグをセット */
					tmps[ii] = 'E';
					ptr++;
					ii++;
				}else{
					/* 既に指数部を読んでいる */
					break;
				}
			}
		}else if(buf[ptr] == '.'){
			/* 小数点らしき物が現れた */
			if(pflag == 1){
				/* 既に読んでいる */
				break;
			}else{
				pflag = 1;
				tmps[ii] = buf[ptr]; /* １文字コピー */
				ptr++;
				ii++;
			}
		}else{
			/* その他の文字が現れた */
			/* 区切りになる文字と同じ動作を行なう */
			if(vflag == 1){
				/* すでに数値らしきものを読んだのでループを出る */
				break;
			}else{
				/* まだ数値らしきものを読んでないのでとりあえず先に進む */
				ptr++;
			}
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}

	tmps[ii] = '\0';

	return(ptr);
}


/* buf[ptr] 以降で整数(らしきもの)をtmpsに切り出す */
static int idc_get_i(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* 数値らしきものが現れたら1 */
	int sflag = 0;  /* 符号が現れたら1 */

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;   /* 文字列の終端に達した */
		}
		if((buf[ptr] == ',')||(buf[ptr] == ' ')||(buf[ptr] == '\t')){
			/* 区切りになる文字が現れた */
			if(vflag == 1){
				/* すでに数値らしきものを読んだのでループを出る */
				break;
			}else{
				/* まだ数値らしきものを読んでないのでとりあえず先に進む */
				ptr++;
			}
		}else if(('0' <= buf[ptr])&&(buf[ptr] <= '9')){
			/* 数字が現れた */
			vflag = 1;
			sflag = 1; /* 符号も現れたことにする */
			tmps[ii] = buf[ptr]; /* １文字コピー */
			ptr++;
			ii++;
		}else if((buf[ptr] == '-')||(buf[ptr] == '+')){
			/* 符号が現れた */
			if(sflag == 1){
				/* 既に符号が現れている */
				break;
			}else{
				/* 符号が現れたのははじめて */
				sflag = 1;
				tmps[ii] = buf[ptr]; /* １文字コピー */
				ptr++;
				ii++;
			}
		}else{
			/* その他の文字が現れた */
			/* 区切りになる文字と同じ動作を行なう */
			if(vflag == 1){
				/* すでに数値らしきものを読んだのでループを出る */
				break;
			}else{
				/* まだ数値らしきものを読んでないのでとりあえず先に進む */
				ptr++;
			}
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}

	tmps[ii] = '\0';

	return(ptr);
}

