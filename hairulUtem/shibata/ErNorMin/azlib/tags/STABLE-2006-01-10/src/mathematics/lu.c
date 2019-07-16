/*********************************************************************
 * lu.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: lu.c 420 2005-08-12 07:02:40Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"

/* Ax = b を解く                     */
/* x と b は同じポインタでも構わない */
/* A の内容は変化しない              */
RC
lu_solve (double **A, int n, const double *b, double *x)
{
	double **LU;
	int *index;
	double d;
	int ii1;

	RC_TRY( allocate2D(n, n, &LU) );
	RC_NULL_CHK( index = mm_alloc(n*sizeof(int)) );

	RC_TRY( copy_matrix(n, n, A, LU) );
	RC_TRY( lu_decomp(LU, n, index, &d) );
	for(ii1=0; ii1<n; ii1++){
		x[ii1] = b[ii1];
	}
	RC_TRY( lu_subst(LU, n, index, x) );

	RC_TRY( free2D(n, n, &LU) );
	RC_TRY( mm_free(index) );

	return(NORMAL_RC);
}


/* 行列 A[n][n] を LU 分解して A に上書き        */
/* index[n] : ピボット選択による行交換の情報     */
/* d : 行交換の回数が偶数なら 1.0, 奇数なら -1.0 */
RC
lu_decomp (double **A, int n, int *index, double *d)
{
	int ii1, ii2, ii3;
	double *scale_factor;
	double max;
	double tmp;
	int max_index;

	if((n <= 0)||(A == NULL)||(index == NULL)||(d == NULL)){
		return(ARG_ERROR_RC);
	}

	scale_factor = (double *)mm_alloc(n * sizeof(double));
	if(scale_factor == NULL){
		return(ALLOC_ERROR_RC);
	}

	/* scale_factor を設定 */
	for(ii1=0; ii1<n; ii1++){
		max = 0.0;
		for(ii2=0; ii2<n; ii2++){
			tmp = fabs(A[ii1][ii2]);
			if(tmp > max) max = tmp;
		}
		scale_factor[ii1] = 1.0/keep_away_zero(max);
	}

	*d = 1.0;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<ii1; ii2++){
			double sum = A[ii2][ii1];

			for(ii3=0; ii3<ii2; ii3++){
				sum -= A[ii2][ii3] * A[ii3][ii1];
			}
			A[ii2][ii1] = sum;
		}

		max = 0.0;
		max_index = ii1;
		for(ii2=ii1; ii2<n; ii2++){
			double sum = A[ii2][ii1];

			for(ii3=0; ii3<ii1; ii3++){
				sum -= A[ii2][ii3] * A[ii3][ii1];
			}
			A[ii2][ii1] = sum;

			tmp = scale_factor[ii2] * fabs(sum);
			if(tmp > max){
				max = tmp;
				max_index = ii2;
			}
		}

		if(max_index != ii1){
			for(ii2=0; ii2<n; ii2++){
				tmp = A[max_index][ii2];
				A[max_index][ii2] = A[ii1][ii2];
				A[ii1][ii2] = tmp;
			}
			*d = -(*d);
			scale_factor[max_index] = scale_factor[ii1];
		}
		index[ii1] = max_index;

		for(ii2=ii1+1; ii2<n; ii2++){
			A[ii2][ii1] /= A[ii1][ii1];
		}
	}

	RC_TRY( mm_free(scale_factor) );

	return(NORMAL_RC);
}


/* 前進代入、後退代入を行なう                                             */
/* A[n][n], index[n] : lu_decomp で得た値                                 */
/* B[n] : 右辺ベクトル -> 解ベクトル                                      */
/* B[] に多くのゼロが含まれる場合(逆マトリックスの計算時)を想定した高速化 */
RC
lu_subst (double **A, int n, int *index, double *B)
{
	int ii1, ii2;
	int nzero_index;
	double tmp;

	if((n <= 0)||(A == NULL)||(index == NULL)||(B == NULL)){
		return(ARG_ERROR_RC);
	}

	/* B[] の並び替え、nzero_index の設定 */
	nzero_index = n;
	for(ii1=0; ii1<n; ii1++){
		if((index[ii1] < 0)||(index[ii1] >= n)) return(ARG_ERROR_RC);

		tmp = B[index[ii1]];
		B[index[ii1]] = B[ii1];
		B[ii1] = tmp;
		if(nzero_index >= n){
			if(!nearly_eq(tmp, 0.0)) nzero_index = ii1;
		}
	}

	/* 前進代入 */
	for(ii1=1; ii1<n; ii1++){
		for(ii2=nzero_index; ii2<ii1; ii2++){
			B[ii1] -= A[ii1][ii2]*B[ii2];
		}
	}

	/* 後退代入 */
	for(ii1=n-1; ii1>=0; ii1--){
		for(ii2=ii1+1; ii2<n; ii2++){
			B[ii1] -= A[ii1][ii2]*B[ii2];
		}
		B[ii1] /= keep_away_zero(A[ii1][ii1]);
	}

	return(NORMAL_RC);
}


