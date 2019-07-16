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


/* buf[ptr] $B0J9_$G6uGrJ8;zEy$G6h@Z$i$l$?J8;zNs$r(Btmps$B$K@Z$j=P$9(B */
static int idc_get_c(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* $BJ8;zNs$,8=$l$?$i(B1 */
	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;  /* $BJ8;zNs(B(buf)$B$N=*C<$KC#$7$?(B */
		}
		if((buf[ptr] == ' ')||(buf[ptr] == '\t')||(buf[ptr] == ',')){
			/* $B6h@Z$j$K$J$kJ8;z$,8=$l$?(B */
			if(vflag == 1){
				/* $B$9$G$KJ8;zNs$rFI$s$@$N$G%k!<%W$r=P$k(B */
				break;
			}else{
				/* $B$^$@J8;zNs$rFI$s$G$J$$$N$G$H$j$"$($:@h$K?J$`(B */
				ptr++;
			}
		}else{
			/* $B0lHLE*$JJ8;z$,8=$l$?(B */
			vflag = 1;
			tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
			ptr++;
			ii++;
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}
	tmps[ii] = '\0';

	return(ptr);
}

/* buf[ptr] $B0J9_$G6uGrJ8;zEy0J30$N#1J8;z$r(Btmps[0]$B$K@Z$j=P$9(B */
static int idc_get_a(const char buf[], char tmps[], int ptr)
{
	int ii = 0;

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;  /* $BJ8;zNs(B(buf)$B$N=*C<$KC#$7$?(B */
		}
		if((buf[ptr] == ' ')||(buf[ptr] == '\t')||(buf[ptr] == ',')){
			/* $B6uGrJ8;z$,8=$l$?(B */
			/* $B@h$K?J$`(B */
			ptr++;
		}else{
			/* $B0lHLE*$JJ8;z$,8=$l$?(B */
			tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
			ptr++;
			ii++;
			break;
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}
	tmps[ii] = '\0';

	return(ptr);
}



/* buf[ptr] $B0J9_$G?tCM(B($B$i$7$-$b$N(B)$B$r(Btmps$B$K@Z$j=P$9(B */
static int idc_get_d(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* $B?tCM$i$7$-$b$N$,8=$l$?$i(B1 */
	int sflag = 0;  /* $BId9f$,8=$l$?$i(B1 */
	int eflag = 0;  /* $B;X?tIt$i$7$-$b$N$,8=$l$?$i(B1 */
	int pflag = 0;  /* $B>.?tE@$i$7$-$b$N$,8=$l$?$i(B1 */

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;   /* $BJ8;zNs$N=*C<$KC#$7$?(B */
		}
		if((buf[ptr] == ',')||(buf[ptr] == ' ')||(buf[ptr] == '\t')){
			/* $B6h@Z$j$K$J$kJ8;z$,8=$l$?(B */
			if(vflag == 1){
				/* $B$9$G$K?tCM$i$7$-$b$N$rFI$s$@$N$G%k!<%W$r=P$k(B */
				break;
			}else{
				/* $B$^$@?tCM$i$7$-$b$N$rFI$s$G$J$$$N$G$H$j$"$($:@h$K?J$`(B */
				ptr++;
			}
		}else if(('0' <= buf[ptr])&&(buf[ptr] <= '9')){
			/* $B?t;z$,8=$l$?(B */
			vflag = 1;
			sflag = 1; /* $BId9f$b8=$l$?$3$H$K$9$k(B */
			tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
			ptr++;
			ii++;
		}else if((buf[ptr] == '-')||(buf[ptr] == '+')){
			/* $BId9f$,8=$l$?(B */
			if(sflag == 1){
				/* $B4{$KId9f$,8=$l$F$$$k(B */
				break;
			}else{
				/* $BId9f$,8=$l$?$N$O$O$8$a$F(B */
				sflag = 1;
				tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
				ptr++;
				ii++;
			}
		}else if((buf[ptr] == 'D')||(buf[ptr] == 'E')
		       ||(buf[ptr] == 'd')||(buf[ptr] == 'e')){
			/* D $B$H$+(B E $B$,8=$l$?(B */
			if(vflag == 0){
				/* $B$^$@?t;z$rFI$s$G$$$J$$$N$GL5;k$7$F<!$K?J$`(B */
				ptr++;
			}else{
				if(eflag == 0){
					eflag = 1; /* $B$3$l$+$i;X?tIt$rFI$`(B */
					sflag = 0; /* $BId9f$N%U%i%0$r$j%;%C%H(B */
					pflag = 1; /* $B>.?tE@$N%U%i%0$r%;%C%H(B */
					tmps[ii] = 'E';
					ptr++;
					ii++;
				}else{
					/* $B4{$K;X?tIt$rFI$s$G$$$k(B */
					break;
				}
			}
		}else if(buf[ptr] == '.'){
			/* $B>.?tE@$i$7$-J*$,8=$l$?(B */
			if(pflag == 1){
				/* $B4{$KFI$s$G$$$k(B */
				break;
			}else{
				pflag = 1;
				tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
				ptr++;
				ii++;
			}
		}else{
			/* $B$=$NB>$NJ8;z$,8=$l$?(B */
			/* $B6h@Z$j$K$J$kJ8;z$HF1$8F0:n$r9T$J$&(B */
			if(vflag == 1){
				/* $B$9$G$K?tCM$i$7$-$b$N$rFI$s$@$N$G%k!<%W$r=P$k(B */
				break;
			}else{
				/* $B$^$@?tCM$i$7$-$b$N$rFI$s$G$J$$$N$G$H$j$"$($:@h$K?J$`(B */
				ptr++;
			}
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}

	tmps[ii] = '\0';

	return(ptr);
}


/* buf[ptr] $B0J9_$G@0?t(B($B$i$7$-$b$N(B)$B$r(Btmps$B$K@Z$j=P$9(B */
static int idc_get_i(const char buf[], char tmps[], int ptr)
{
	int ii = 0;
	int vflag = 0;  /* $B?tCM$i$7$-$b$N$,8=$l$?$i(B1 */
	int sflag = 0;  /* $BId9f$,8=$l$?$i(B1 */

	while(1){
		if((buf[ptr] == '\0')||(buf[ptr] == '\n')){
			break;   /* $BJ8;zNs$N=*C<$KC#$7$?(B */
		}
		if((buf[ptr] == ',')||(buf[ptr] == ' ')||(buf[ptr] == '\t')){
			/* $B6h@Z$j$K$J$kJ8;z$,8=$l$?(B */
			if(vflag == 1){
				/* $B$9$G$K?tCM$i$7$-$b$N$rFI$s$@$N$G%k!<%W$r=P$k(B */
				break;
			}else{
				/* $B$^$@?tCM$i$7$-$b$N$rFI$s$G$J$$$N$G$H$j$"$($:@h$K?J$`(B */
				ptr++;
			}
		}else if(('0' <= buf[ptr])&&(buf[ptr] <= '9')){
			/* $B?t;z$,8=$l$?(B */
			vflag = 1;
			sflag = 1; /* $BId9f$b8=$l$?$3$H$K$9$k(B */
			tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
			ptr++;
			ii++;
		}else if((buf[ptr] == '-')||(buf[ptr] == '+')){
			/* $BId9f$,8=$l$?(B */
			if(sflag == 1){
				/* $B4{$KId9f$,8=$l$F$$$k(B */
				break;
			}else{
				/* $BId9f$,8=$l$?$N$O$O$8$a$F(B */
				sflag = 1;
				tmps[ii] = buf[ptr]; /* $B#1J8;z%3%T!<(B */
				ptr++;
				ii++;
			}
		}else{
			/* $B$=$NB>$NJ8;z$,8=$l$?(B */
			/* $B6h@Z$j$K$J$kJ8;z$HF1$8F0:n$r9T$J$&(B */
			if(vflag == 1){
				/* $B$9$G$K?tCM$i$7$-$b$N$rFI$s$@$N$G%k!<%W$r=P$k(B */
				break;
			}else{
				/* $B$^$@?tCM$i$7$-$b$N$rFI$s$G$J$$$N$G$H$j$"$($:@h$K?J$`(B */
				ptr++;
			}
		}
		if(ii >= (BUFFER_SIZE - 1))break;
	}

	tmps[ii] = '\0';

	return(ptr);
}

