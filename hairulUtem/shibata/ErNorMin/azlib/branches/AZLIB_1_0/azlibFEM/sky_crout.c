/*********************************************************************
 * sky_crout.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sky_crout.c,v 1.6 2003/07/23 04:06:32 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"

/*
0 2   A
1 3 5 B
  4 6 C
7 8 9 D
*/

#define NE_ZERO    (1.0e-50)  /* 分母の絶対値がこれより小さければエラー */

/* 配列の確保 */
/* index1[] を参照し、index2[],array[] を確保 */
/* 確保した配列は適当に初期化される */
/* index2[] の長さは n */
/* array_size が NULL でなければ、array[] の長さが array_size に代入される */
RC s_crout_alloc(long n, const long *index1, long **index2,
                 double **array, long *array_size)
{
	long ii1;
	long sum;
	
	if(n < 1L)return(ARG_ERROR_RC);

	/* index2[] を確保 */
	*index2 = (long *)malloc(n*sizeof(long));
	if((*index2) == (long *)NULL)return(ALLOC_ERROR_RC);

	/* array[] の要素数をカウント */
	/* index2[] の初期化 */
	sum = 0L;
	for(ii1=0L;ii1<n;ii1++){
		sum += 2L*(ii1 - index1[ii1]) + 1L;
		(*index2)[ii1] = sum - 1L;
	}
	if(array_size != NULL) *array_size = sum;
	
	/* array[] を確保 */
	*array = (double *)malloc(sum*sizeof(double));
	if((*array) == (double *)NULL)return(ALLOC_ERROR_RC);

	for(ii1=0L; ii1<sum; ii1++){
		(*array)[ii1] = 0.0;
	}
	
	return(NORMAL_RC);
}


/* 配列の解放 */
RC s_crout_free(long *index2, double *array)
{
	if(((void *)array == NULL)||((void *)index2 == NULL)){
		return(ARG_ERROR_RC);
	}
	
	free((void *)array);
	free((void *)index2);
	
	return(NORMAL_RC);
}


double s_crout_get(long n, const long *index1,
                           const long *index2, const double *array,
						   int row, int col, RC *rc)
{
	long ix;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)){
		(*rc) = ARG_ERROR_RC;
		return(0.0);
	}

	if(row > col){
		if(col < index1[row]){
			(*rc) = SEEK_ERROR_RC;
			return(0.0);
		}
		ix = index2[row] -2L*row + index1[row] + col;
	}else{
		if(row < index1[col]){
			(*rc) = SEEK_ERROR_RC;
			return(0.0);
		}
		ix = index2[col] - col + row;
	}

	(*rc) = NORMAL_RC;
	return(array[ix]);
}


double *s_crout_ptr(long n, const long *index1,
                            const long *index2, double *array,
                            int row, int col, RC *rc)
{
	long ix;
	static double dummy;

	dummy = 0.0;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)){
		(*rc) = ARG_ERROR_RC;
		return(&dummy);
	}

	if(row > col){
		if(col < index1[row]){
			(*rc) = SEEK_ERROR_RC;
			return(&dummy);
		}
		ix = index2[row] -2L*row + index1[row] + col;
	}else{
		if(row < index1[col]){
			(*rc) = SEEK_ERROR_RC;
			return(&dummy);
		}
		ix = index2[col] - col + row;
	}

	(*rc) = NORMAL_RC;
	return( &(array[ix]) );
}


/* 配列の内容を表示 (debug用) */
void s_crout_print(long n, const long *index1,
                           const long *index2, const double *array)
{
    long ii1,ii2;
	RC rc;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)){
		fprintf(stderr, "???\n");
		return;
	}

    for(ii1=0L;ii1<n;ii1++){
		printf("%ld  %ld\n",index1[ii1],index2[ii1]);
        for(ii2=0L;ii2<n;ii2++){
            printf("%10.3f ",s_crout_get(n,index1,index2,array,ii1,ii2,&rc));
        }
        printf("\n");
    }
}


/* m 番目の解が常に v となるよう、array[],vector[] を修正 */
RC s_crout_mod(long n, long *index1, const long *index2,
                          double *array, double *vector, long m, double v)
{
	long ii1;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
	           ||((void *)vector == NULL)
	           ||(n <= m)) return(ARG_ERROR_RC);

	for(ii1=index1[m];ii1<m;ii1++){
		vector[ii1] -= v * array[index2[m] - m + ii1];
		array[index2[m] - m + ii1] = 0.0;
		array[index2[m] -2L*m + index1[m] + ii1] = 0.0;
	}

	array[index2[m]] = 1.0;
	index1[m] = m;

	for(ii1=m+1;ii1<n;ii1++){
		if(index1[ii1] <= m){
			vector[ii1] -= v * array[index2[ii1] -2L*ii1 + index1[ii1] + m];
			array[index2[ii1] - ii1 + m] = 0.0;
			array[index2[ii1] - 2L*ii1 + index1[ii1] + m] = 0.0;
			/* Do not do it !! */
			/* if((index1[ii1] == m)&&(index1[ii1] > ii1)) index1[ii1]--; */
		}
	}
	vector[m] = v;

	return(NORMAL_RC);
}


/* クラウト法による分解 */
RC s_crout_decomp(long n, const long *index1, const long *index2,
                                              double *array, int p_flag)
{
	long ii1, ii2, ii3;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)) return(ARG_ERROR_RC);
	
	for(ii1=0L;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<=ii1;ii2++){
			if(index1[ii1] > index1[ii2]){
				ii3=index1[ii1];
			}else{
				ii3=index1[ii2];
			}
			for( ;ii3<ii2;ii3++){
				array[index2[ii1] - ii1 + ii2]
				  -= array[index2[ii1] - ii1 + ii3]
				   * array[index2[ii2] - 2L*ii2 + index1[ii2] + ii3];
			}
		}

		if(fabs(array[index2[ii1]]) < NE_ZERO){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}

		for(ii2=ii1+1;ii2<n;ii2++){
			if(index1[ii2] > ii1) continue;
			if(index1[ii1] > index1[ii2]){
				ii3=index1[ii1];
			}else{
				ii3=index1[ii2];
			}
			for( ;ii3<ii1;ii3++){
				array[index2[ii2] - 2L*ii2 + index1[ii2] + ii1]
				  -= array[index2[ii2] - 2L*ii2 + index1[ii2] + ii3]
				   * array[index2[ii1] - ii1 + ii3];
			}

			array[index2[ii2] - 2L*ii2 + index1[ii2] + ii1]
			                          /= array[index2[ii1]];
		}
	}

	return(NORMAL_RC);
}


/* 前進代入、後退代入 */
RC s_crout_subst(long n, const long *index1, const long *index2,
                         const double *array, double *vector, int p_flag)
{
	int ii1, ii2;


	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
			   ||((void *)vector == NULL)) return(ARG_ERROR_RC);

	/* 前進代入 */
	for(ii1=0;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii1] -= vector[ii2]
			             * array[index2[ii1] - 2L*ii1 + index1[ii1] + ii2];
		}
	}

	/* 後退代入 */
	for(ii1=n-1;ii1>=0;ii1--){
		vector[ii1] /= array[index2[ii1]];

		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii2] -= vector[ii1] * array[index2[ii1] - ii1 + ii2];
		}
	}

	return(NORMAL_RC);
}



