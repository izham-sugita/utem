/*********************************************************************
 * random.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Masanobu SUYAMA> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: random.c 1075 2016-10-31 09:48:59Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "rc.h"
#include "math_utl.h"

/* NUMERICAL RECIPES in C より */

/* 乱数用パラメータ */
#define D_RANDOM_IM1 2147483563
#define D_RANDOM_IM2 2147483399
#define D_RANDOM_AM (1.0/D_RANDOM_IM1) /* 乱数列の最小値 */
#define D_RANDOM_IMM1 (D_RANDOM_IM1-1)
#define D_RANDOM_IA1 40014
#define D_RANDOM_IA2 40692
#define D_RANDOM_IQ1 53668
#define D_RANDOM_IQ2 52774
#define D_RANDOM_IR1 12211
#define D_RANDOM_IR2 3791
#define D_RANDOM_NTAB 32
#define D_RANDOM_NDIV (1+D_RANDOM_IMM1/D_RANDOM_NTAB)
#define D_RANDOM_EPS 1.2e-7
#define D_RANDOM_RNMX (1.0-D_RANDOM_EPS) /* 乱数列の最大値 */

static long idum1 = 1L;
static long idum2 = 123456789L;
static long iy = 0L;
static long iv[D_RANDOM_NTAB];

/* 乱数列を idum で初期化する  */
void
init_d_random (long idum)
{
	int j;
	long k;

	if(idum < 1){
		idum1 = 1;
	}else{
		idum1 = idum;
	}
	idum2 = 123456789;
	for(j=D_RANDOM_NTAB+7;j>=0;j--){
		k = idum1/D_RANDOM_IQ1;
		idum1 = D_RANDOM_IA1*(idum1 - k*D_RANDOM_IQ1) - k*D_RANDOM_IR1;
		if(idum1 < 0) idum1 += D_RANDOM_IM1;
		if(j < D_RANDOM_NTAB)iv[j] = idum1;
	}
	iy = iv[0];
}


/* 0 ~ 1 の乱数列を発生させる */
double
d_random (void)
{
	int j;
	long k;
	double temp;

	if(iy <= 0){
		init_d_random(1L);
	}

	k = (idum1)/D_RANDOM_IQ1;
	idum1 = D_RANDOM_IA1*(idum1 - k*D_RANDOM_IQ1) - k*D_RANDOM_IR1;
	if(idum1 < 0) idum1 += D_RANDOM_IM1;

	k = idum2/D_RANDOM_IQ2;
	idum2 = D_RANDOM_IA2*(idum2 - k*D_RANDOM_IQ2) - k*D_RANDOM_IR2;
	if(idum2 < 0)idum2 += D_RANDOM_IM2;

	j = iy/D_RANDOM_NDIV;
	iy = iv[j] - idum2;
	iv[j] = idum1;
	if(iy < 1)iy += D_RANDOM_IMM1;

	if((temp=D_RANDOM_AM*iy) > D_RANDOM_RNMX){
	 	return D_RANDOM_RNMX;
	}
	return temp;
}


/* 現在の時間を乱数の種として返却          */
/* 62*60*24*31*12*(2005-1900) = 3487276800 */
/* これぐらいの大きな数字になる            */
long
random_seed_time (void)
{
	time_t tt;
	struct tm *tm_ptr;
    
	tt = time(NULL);
	if(tt == (time_t)-1){
		warning_printf(0, "time() error.\n");
		return(0);
	}

	tm_ptr = localtime(&tt);
	if(!tm_ptr){
		warning_printf(0, "localtime() error.\n");
		return(0);
	}

	return( (long)(tm_ptr->tm_sec+1)*(tm_ptr->tm_min+1)*(tm_ptr->tm_hour+1)
	        *(tm_ptr->tm_mday)*(tm_ptr->tm_mon+1)*(tm_ptr->tm_year));

}


/* 配列 array の n 個の要素にランダムな実数を代入する．
 */
RC
fill_random_array (int n, double array[], int seed)
{
	int ii1;

	init_d_random((long)seed);
	for(ii1=0; ii1<n; ii1++){
		array[ii1] = d_random();
	}

	return(NORMAL_RC);
}


