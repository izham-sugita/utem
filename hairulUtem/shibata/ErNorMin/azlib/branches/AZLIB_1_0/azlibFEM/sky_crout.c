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

#define NE_ZERO    (1.0e-50)  /* $BJ,Jl$N@dBPCM$,$3$l$h$j>.$5$1$l$P%(%i!<(B */

/* $BG[Ns$N3NJ](B */
/* index1[] $B$r;2>H$7!"(Bindex2[],array[] $B$r3NJ](B */
/* $B3NJ]$7$?G[Ns$OE,Ev$K=i4|2=$5$l$k(B */
/* index2[] $B$ND9$5$O(B n */
/* array_size $B$,(B NULL $B$G$J$1$l$P!"(Barray[] $B$ND9$5$,(B array_size $B$KBeF~$5$l$k(B */
RC s_crout_alloc(long n, const long *index1, long **index2,
                 double **array, long *array_size)
{
	long ii1;
	long sum;
	
	if(n < 1L)return(ARG_ERROR_RC);

	/* index2[] $B$r3NJ](B */
	*index2 = (long *)malloc(n*sizeof(long));
	if((*index2) == (long *)NULL)return(ALLOC_ERROR_RC);

	/* array[] $B$NMWAG?t$r%+%&%s%H(B */
	/* index2[] $B$N=i4|2=(B */
	sum = 0L;
	for(ii1=0L;ii1<n;ii1++){
		sum += 2L*(ii1 - index1[ii1]) + 1L;
		(*index2)[ii1] = sum - 1L;
	}
	if(array_size != NULL) *array_size = sum;
	
	/* array[] $B$r3NJ](B */
	*array = (double *)malloc(sum*sizeof(double));
	if((*array) == (double *)NULL)return(ALLOC_ERROR_RC);

	for(ii1=0L; ii1<sum; ii1++){
		(*array)[ii1] = 0.0;
	}
	
	return(NORMAL_RC);
}


/* $BG[Ns$N2rJ|(B */
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


/* $BG[Ns$NFbMF$rI=<((B (debug$BMQ(B) */
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


/* m $BHVL\$N2r$,>o$K(B v $B$H$J$k$h$&!"(Barray[],vector[] $B$r=$@5(B */
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


/* $B%/%i%&%HK!$K$h$kJ,2r(B */
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


/* $BA0?JBeF~!"8eB`BeF~(B */
RC s_crout_subst(long n, const long *index1, const long *index2,
                         const double *array, double *vector, int p_flag)
{
	int ii1, ii2;


	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
			   ||((void *)vector == NULL)) return(ARG_ERROR_RC);

	/* $BA0?JBeF~(B */
	for(ii1=0;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii1] -= vector[ii2]
			             * array[index2[ii1] - 2L*ii1 + index1[ii1] + ii2];
		}
	}

	/* $B8eB`BeF~(B */
	for(ii1=n-1;ii1>=0;ii1--){
		vector[ii1] /= array[index2[ii1]];

		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii2] -= vector[ii1] * array[index2[ii1] - ii1 + ii2];
		}
	}

	return(NORMAL_RC);
}



