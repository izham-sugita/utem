/*********************************************************************
 * sky_cholesky.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: sky_cholesky.c,v 1.8 2003/07/23 04:06:32 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"

/*
0 1   5
1 2 3 6
  3 4 7
5 6 7 8
*/

#define NE_ZERO    (1.0e-50)  /* $BJ,Jl$N@dBPCM$,$3$l$h$j>.$5$1$l$P%(%i!<(B */

/* $BG[Ns$N3NJ](B */
/* index1[] $B$r;2>H$7!"(Bindex2[],array[] $B$r3NJ](B */
/* $B3NJ]$7$?G[Ns$OE,Ev$K=i4|2=$5$l$k(B */
/* index2[] $B$ND9$5$O(B n */
/* array_size $B$,(B NULL $B$G$J$1$l$P!"(Barray[] $B$ND9$5$,(B array_size $B$KBeF~$5$l$k(B */
RC s_cholesky_alloc(long n, const long *index1, long **index2,
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
		sum += ii1 - index1[ii1] + 1L;
		(*index2)[ii1] = sum - 1L;
	}
	if(array_size != NULL) *array_size = sum;
#ifdef DEBUG
	fprintf(stderr, "s_cholesky.c: allocated size: %ld\n", sum);
#endif  /* DEBUG */
	
	/* array[] $B$r3NJ](B */
	*array = (double *)malloc(sum*sizeof(double));
	if((*array) == (double *)NULL)return(ALLOC_ERROR_RC);

	for(ii1=0L; ii1<sum; ii1++){
		(*array)[ii1] = 0.0;
	}
	
	return(NORMAL_RC);
}


/* $BG[Ns$N2rJ|(B */
RC s_cholesky_free(long *index2, double *array)
{
	if(((void *)array == NULL)||((void *)index2 == NULL)){
		return(ARG_ERROR_RC);
	}
	
	free((void *)array);
	free((void *)index2);
	
	return(NORMAL_RC);
}


double s_cholesky_get(long n, const long *index1, const long *index2,
                      const double *array, int row, int col, RC *rc)
{
	long ix;
	int tmp;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)){
		(*rc) = ARG_ERROR_RC;
		return(0.0);
	}

	if(row > col){
		tmp = col;
		col = row;
		row = tmp;
	}

	if(row < index1[col]){
		(*rc) = SEEK_ERROR_RC;
		return(0.0);
	}
	ix = index2[col] - col + row;

	(*rc) = NORMAL_RC;
	return(array[ix]);
}


double *s_cholesky_ptr(long n, const long *index1, const long *index2,
                       double *array, int row, int col, RC *rc)
{
	long ix;
	int tmp;
	static double dummy;

	dummy = 0.0;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)){
		(*rc) = ARG_ERROR_RC;
		return(&dummy);
	}

	if(row > col){
		tmp = col;
		col = row;
		row = tmp;
	}

	if(row < index1[col]){
		(*rc) = SEEK_ERROR_RC;
		return(&dummy);
	}
	ix = index2[col] - col + row;

	(*rc) = NORMAL_RC;
	return( &(array[ix]) );
}


/* $BG[Ns$NFbMF$rI=<((B (debug$BMQ(B) */
void s_cholesky_print(long n, const long *index1,
                              const long *index2, const double *array)
{
    long ii1,ii2;

    for(ii1=0L;ii1<n;ii1++){
		printf("%ld  %ld\n", index1[ii1], index2[ii1]);
        for(ii2=index2[ii1]+index1[ii1]-ii1; ii2<=index2[ii1]; ii2++){
            printf("%f  ", array[ii2]);
        }
        printf("\n");
    }
}


/* m $BHVL\$N2r$,>o$K(B v $B$H$J$k$h$&!"(Barray[],vector[] $B$r=$@5(B */
RC s_cholesky_mod(long n, long *index1, const long *index2,
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
	}
	array[index2[m]] = 1.0;
	index1[m] = m;

	for(ii1=m+1;ii1<n;ii1++){
		if(index1[ii1] <= m){
			vector[ii1] -= v * array[index2[ii1] - ii1 + m];
			array[index2[ii1] - ii1 + m] = 0.0;
			if(index1[ii1] == m) index1[ii1]++;
		}
	}
	vector[m] = v;

	return(NORMAL_RC);
}


/* m $BHVL\$N2r$,>o$K(B v $B$H$J$k$h$&!"(Barray[],vectors[vec_size][] $B$r=$@5(B */
RC s_cholesky_mod_n(long n, long *index1, const long *index2,
                    double *array, double **vectors, long m, double v,
                    int vec_size)
{
	long ii1, ii2;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
	           ||((void *)vectors == NULL)
	           ||(n <= m)) return(ARG_ERROR_RC);

	for(ii1=index1[m];ii1<m;ii1++){
		for(ii2=0; ii2<vec_size; ii2++){
			vectors[ii2][ii1] -= v * array[index2[m] - m + ii1];
		}
		array[index2[m] - m + ii1] = 0.0;
	}
	array[index2[m]] = 1.0;
	index1[m] = m;

	for(ii1=m+1;ii1<n;ii1++){
		if(index1[ii1] <= m){
			for(ii2=0; ii2<vec_size; ii2++){
				vectors[ii2][ii1] -= v * array[index2[ii1] - ii1 + m];
			}
			array[index2[ii1] - ii1 + m] = 0.0;
			if(index1[ii1] == m) index1[ii1]++;
		}
	}

	for(ii2=0; ii2<vec_size; ii2++){
		vectors[ii2][m] = v;
	}

	return(NORMAL_RC);
}


#if 0
/* $B2~D{%3%l%9%-!<K!$K$h$kJ,2r(B */
RC s_cholesky_decomp(long n, const long *index1, const long *index2,
                                               double *array, int p_flag)
{
	long ii1, ii2, ii3;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)) return(ARG_ERROR_RC);

	/*s_cholesky_print(n,index1,index2,array);*/

	for(ii1=0L;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<ii1;ii2++){

			ii3 = index1[ii1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];

			for(;ii3<ii2;ii3++){
				array[index2[ii1] - ii1 + ii2]
				   -= array[index2[ii1] - ii1 + ii3]
					* array[index2[ii2] - ii2 + ii3]
				    * array[index2[ii3]];
			}
			array[index2[ii1] - ii1 + ii2] /= array[index2[ii2]];
			array[index2[ii1]] 
			     -= array[index2[ii1] - ii1 + ii2]
				  * array[index2[ii1] - ii1 + ii2]
				  * array[index2[ii2]];
		}
		if(fabs(array[index2[ii1]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}
	}

	return(NORMAL_RC);
}
#endif /* 0 */


/* $B2~D{%3%l%9%-!<K!$K$h$kJ,2r(B($B9bB.HG(B) */
RC s_cholesky_decomp(long n, const long *index1, const long *index2,
                                               double *array, int p_flag)
{
	long ii1, ii2, ii3;
	double *mul;    /* $B>h;;2s?tDc8:MQ(B */

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)) return(ARG_ERROR_RC);

	/* s_cholesky_print(n,index1,index2,array); */

	/* mul[] $B$r3NJ](B */
	mul = (double *)malloc(n*sizeof(double));
	if(mul == (double *)NULL)return(ALLOC_ERROR_RC);

	for(ii1=0L;ii1<n;ii1++){
		if(p_flag){
			fprintf(stderr, " [F%6ld/%6ld]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
			        ii1, n);
		}
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			double tmp = array[index2[ii1] - ii1 + ii2];
			long shift_ii2 = index2[ii2] - ii2;

			ii3 = index1[ii1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];

			for(;ii3<ii2;ii3++){
				tmp -= mul[ii3] * array[shift_ii2 + ii3];
			}
			mul[ii2] = tmp;
			tmp /= array[index2[ii2]];
			array[index2[ii1]] -= tmp * tmp * array[index2[ii2]];

			array[index2[ii1] - ii1 + ii2] = tmp;
		}
		if(fabs(array[index2[ii1]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}
	}

	if(p_flag) fprintf(stderr, " [F%6ld/%6ld]", n, n);

	free((void *)mul);

	return(NORMAL_RC);
}


/* $B2~D{%3%l%9%-!<K!$K$h$kJ,2r(B($B9bB.HG(B) */
RC s_cholesky_decomp3(long n, const long *index1, const long *index2,
                                               double *array, int p_flag)
{
	long ii1, ii2, ii3;
	long jj1;
	double *mul;    /* $B>h;;2s?tDc8:MQ(B */

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)) return(ARG_ERROR_RC);

	/* s_cholesky_print(n,index1,index2,array); */

	/* mul[] $B$r3NJ](B */
	mul = (double *)malloc(3*n*sizeof(double));
	if(mul == (double *)NULL)return(ALLOC_ERROR_RC);

	for(ii1=0L;ii1<(n-2);ii1+=3){
		int max_index;

		if(p_flag&&(ii1%120L == 0)){
			fprintf(stderr, " [F%6ld/%6ld]\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b",
			        ii1, n);
		}

		max_index = index1[ii1];
		if(index1[ii1+1] > max_index) max_index = index1[ii1+1];
		if(index1[ii1+2] > max_index) max_index = index1[ii1+2];
		if(max_index > ii1) max_index = ii1;

		for(jj1=0; jj1<3; jj1++){
			for(ii2=index1[ii1+jj1]; ii2<max_index; ii2++){
				double tmp = array[index2[ii1+jj1] - (ii1+jj1) + ii2];

				ii3 = index1[ii1+jj1];
				if(ii3 < index1[ii2])ii3 = index1[ii2];

				for(;ii3<ii2;ii3++){
					tmp -= mul[ii3*3 + jj1] * array[index2[ii2] - ii2 + ii3];
				}
				mul[ii2*3 + jj1] = tmp;
				tmp /= array[index2[ii2]];
				array[index2[ii1+jj1]] -= tmp * tmp * array[index2[ii2]];

				array[index2[ii1+jj1] - (ii1+jj1) + ii2] = tmp;
			}
		}

		for(ii2=max_index;ii2<ii1;ii2++){
			double tmp0 = array[index2[ii1] - ii1 + ii2];
			double tmp1 = array[index2[ii1+1] - (ii1+1) + ii2];
			double tmp2 = array[index2[ii1+2] - (ii1+2) + ii2];
			long shift_ii2 = index2[ii2] - ii2;

			ii3 = index1[ii1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];
			for(;ii3<max_index;ii3++){
				tmp0 -= mul[ii3*3] * array[shift_ii2 + ii3];
			}

			ii3 = index1[ii1 + 1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];
			for(;ii3<max_index;ii3++){
				tmp1 -= mul[ii3*3 + 1] * array[shift_ii2 + ii3];
			}

			ii3 = index1[ii1 + 2];
			if(ii3 < index1[ii2])ii3 = index1[ii2];
			for(;ii3<max_index;ii3++){
				tmp2 -= mul[ii3*3 + 2] * array[shift_ii2 + ii3];
			}

			ii3 = max_index;
			if(ii3 < index1[ii2])ii3 = index1[ii2];
			for(;ii3<ii2;ii3++){
				double v = array[shift_ii2 + ii3];
				tmp0 -= mul[ii3*3] * v;
				tmp1 -= mul[ii3*3+1] * v;
				tmp2 -= mul[ii3*3+2] * v;
			}
			mul[ii2*3] = tmp0;
			mul[ii2*3+1] = tmp1;
			mul[ii2*3+2] = tmp2;
			tmp0 /= array[index2[ii2]];
			tmp1 /= array[index2[ii2]];
			tmp2 /= array[index2[ii2]];
			array[index2[ii1]] -= tmp0 * tmp0 * array[index2[ii2]];
			array[index2[ii1+1]] -= tmp1 * tmp1 * array[index2[ii2]];
			array[index2[ii1+2]] -= tmp2 * tmp2 * array[index2[ii2]];

			array[index2[ii1] - ii1 + ii2] = tmp0;
			array[index2[ii1+1] - (ii1+1) + ii2] = tmp1;
			array[index2[ii1+2] - (ii1+2) + ii2] = tmp2;
		}
		if(fabs(array[index2[ii1]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}

		ii2 = ii1;
		if(index1[ii1+1] > ii2) ii2 = index1[ii1+1];
		for(;ii2<(ii1+1);ii2++){
			double tmp = array[index2[ii1+1] - (ii1+1) + ii2];

			ii3 = index1[ii1+1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];

			for(;ii3<ii2;ii3++){
				tmp -= mul[ii3*3 + 1] * array[index2[ii2] - ii2 + ii3];
			}
			mul[ii2*3 +1] = tmp;
			tmp /= array[index2[ii2]];
			array[index2[ii1+1]] -= tmp * tmp * array[index2[ii2]];

			array[index2[ii1+1] - (ii1+1) + ii2] = tmp;
		}
		if(fabs(array[index2[ii1+1]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}

		ii2 = ii1;
		if(index1[ii1+2] > ii2) ii2 = index1[ii1+2];
		for(;ii2<(ii1+2);ii2++){
			double tmp = array[index2[ii1+2] - (ii1+2) + ii2];

			ii3 = index1[ii1+2];
			if(ii3 < index1[ii2])ii3 = index1[ii2];

			for(;ii3<ii2;ii3++){
				tmp -= mul[ii3*3 + 2] * array[index2[ii2] - ii2 + ii3];
			}
			mul[ii2*3 + 2] = tmp;
			tmp /= array[index2[ii2]];
			array[index2[ii1+2]] -= tmp * tmp * array[index2[ii2]];

			array[index2[ii1+2] - (ii1+2) + ii2] = tmp;
		}
		if(fabs(array[index2[ii1+2]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}
	}

	for(;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			double tmp = array[index2[ii1] - ii1 + ii2];

			ii3 = index1[ii1];
			if(ii3 < index1[ii2])ii3 = index1[ii2];

			for(;ii3<ii2;ii3++){
				tmp -= mul[ii3*3] * array[index2[ii2] - ii2 + ii3];
			}
			mul[ii2*3] = tmp;
			tmp /= array[index2[ii2]];
			array[index2[ii1]] -= tmp * tmp * array[index2[ii2]];

			array[index2[ii1] - ii1 + ii2] = tmp;
		}
		if(fabs(array[index2[ii1]]) < ABS_TOL){
			if(p_flag)fprintf(stderr,"pivoting error\n");
			return(CAL_ERROR_RC);
		}
	}

	if(p_flag) fprintf(stderr, " [F%6ld/%6ld]", n, n);

	free((void *)mul);

	return(NORMAL_RC);
}


/* $B2~D{%3%l%9%-!<K!$K$h$jF@$i$l$?2<;03Q9TNs!"BP3Q9TNs$rMQ$$$F(B */
/* $BA0?JBeF~!"8eB`BeF~$r9T$J$&(B */
RC s_cholesky_subst(long n, const long *index1, const long *index2,
                            const double *array, double *vector, int p_flag)
{
	long ii1, ii2;


	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
			   ||((void *)vector == NULL)) return(ARG_ERROR_RC);

	/* $BA0?JBeF~(B */
	if(p_flag) fprintf(stderr, "[B]");
	for(ii1=0L;ii1<n;ii1++){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii1] -= vector[ii2] * array[index2[ii1] - ii1 + ii2];
		}
	}

	/* $BBP3Q9TNs(B */
	for(ii1=0L;ii1<n;ii1++){
		if(fabs(array[index2[ii1]]) < ABS_TOL){
			return(CAL_ERROR_RC);
		}
		vector[ii1] /= array[index2[ii1]];
	}

	/* $B8eB`BeF~(B */
	for(ii1=n-1L;ii1>=0;ii1--){
		for(ii2=index1[ii1];ii2<ii1;ii2++){
			vector[ii2] -= vector[ii1] * array[index2[ii1] - ii1 + ii2];
		}
	}
	if(p_flag) fprintf(stderr, "...finished.\n");

	return(NORMAL_RC);
}


/* $B9TNs$H(B vector $B$r3]$1$F(B answer $B$KBeF~(B */
RC s_cholesky_mul_vect(long n, const long *index1, const long *index2,
                       const double *array, const double *vector,
                       double *answer)
{
	int ii1;
	int ii2;

	if((n < 1L)||((void *)array == NULL)
	           ||((void *)index1 == NULL)
	           ||((void *)index2 == NULL)
			   ||((void *)vector == NULL)
	           ||((void *)answer == NULL)) return(ARG_ERROR_RC);

	for(ii1=0; ii1<n; ii1++){
		answer[ii1] = 0.0;
		for(ii2=index1[ii1]; ii2<ii1; ii2++){
			answer[ii1] += vector[ii2] * array[index2[ii1] - ii1 + ii2];
			answer[ii2] += vector[ii1] * array[index2[ii1] - ii1 + ii2];
		}
		answer[ii1] += vector[ii1] * array[index2[ii2]];
	}

	return(NORMAL_RC);
}


