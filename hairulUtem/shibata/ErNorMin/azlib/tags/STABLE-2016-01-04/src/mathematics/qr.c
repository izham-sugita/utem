/*********************************************************************
 * qr.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: qr.c 842 2006-12-04 06:21:20Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"


RC
qr_solve (double **A, int n, const double *b, double *x)
{
	double **Q;
	double **R;

	RC_TRY( allocate2D(n, n, &Q) );
	RC_TRY( allocate2D(n, n, &R) );

	RC_TRY( qr_factor(n, A, Q, R) );
	RC_TRY( qr_subst(n, Q, R, b, x) );

	RC_TRY( free2D(n, n, &Q) );
	RC_TRY( free2D(n, n, &R) );

	return(NORMAL_RC);
}


/* A[n][n] を直交行列 Q[n][n] と 上(右)三角R[n][n] に分解 */
RC
qr_factor (int n, double **A, double **Q, double **R)
{
	int ii1, ii2, ii3, ii4;
	double tmp_product, max_product;
	double min_norm;
	int *zero_flags;
	double *vect_sum;
	int max_index;

	RC_NULL_CHK( zero_flags = (int *)mm_alloc(sizeof(int)*n) );
	for(ii1=0; ii1<n; ii1++) zero_flags[ii1] = 0;
	RC_NULL_CHK( vect_sum = (double *)mm_alloc(sizeof(double)*n) );

	init_matrix(n, n, R);

	/* 列同士の内積を高速化するため A^t => Q */
	/* 以後，A の列参照は Q の行参照で代用 */
	transpose_matrix(n, n, A, Q);

	min_norm = ABS_TOL;
	for(ii1=0; ii1<n; ii1++){
		tmp_product = sqrt( inner_product(n, Q[ii1], Q[ii1]) );
		if(tmp_product > min_norm) min_norm = tmp_product;
	}
	min_norm *= (10e+3)*REL_TOL;

	for(ii1=0; ii1<n; ii1++){
		/* G-S 直交化 */
		for(ii2=0; ii2<ii1; ii2++){
			R[ii2][ii1] = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				Q[ii1][ii3] -= R[ii2][ii1]*Q[ii2][ii3];
			}
		}

		/* G-S 直交化を再度行って精度向上 */
		for(ii2=0; ii2<n; ii2++) vect_sum[ii2] = 0.0;  /* 桁落ち防止 */
		for(ii2=0; ii2<ii1; ii2++){
			tmp_product = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				vect_sum[ii3] -= tmp_product*Q[ii2][ii3];  /* 桁落ち防止 */
			}
			R[ii2][ii1] += tmp_product;
		}
		for(ii2=0; ii2<n; ii2++) Q[ii1][ii2] += vect_sum[ii2]; /* 桁落ち防止 */

		R[ii1][ii1] = sqrt( inner_product(n, Q[ii1], Q[ii1]) );
		if(R[ii1][ii1] > min_norm){
			for(ii2=0; ii2<n; ii2++){
				Q[ii1][ii2] /= R[ii1][ii1];
			}
		}else{
			zero_flags[ii1] = 1;
			for(ii2=0; ii2<n; ii2++){
				Q[ii1][ii2] = 0.0;
			}
		}
	}

	/* zero_flags の立っている R の対角項をゼロにする */
	for(ii1=0; ii1<n; ii1++){
		if(zero_flags[ii1] == 0) continue;
		R[ii1][ii1] = 0.0;
	}

	/* zero_flags の立っている Q の行を適当なベクトルで埋める */
	for(ii1=0; ii1<n; ii1++){
		if(zero_flags[ii1] == 0) continue;
		max_product = -1.0;
		max_index = 0;
		for(ii2=0; ii2<n; ii2++){
			/* ベクトルの候補を設定 */
			for(ii3=0; ii3<n; ii3++){
				Q[ii1][ii3] = 0.0;
			}
			Q[ii1][ii2] = 1.0;

			/* G-S 直交化 */
			for(ii3=0; ii3<n; ii3++){
				if(zero_flags[ii3]) continue;
				tmp_product = inner_product(n, Q[ii1], Q[ii3]);
				for(ii4=0; ii4<n; ii4++){
					Q[ii1][ii4] -= tmp_product*Q[ii3][ii4];
				}
			}
			tmp_product = sqrt( inner_product(n, Q[ii1], Q[ii1]) );
			if(tmp_product > max_product){
				max_index = ii2;
				max_product = tmp_product;
			}
		}

		/* ベクトル設定 */
		for(ii2=0; ii2<n; ii2++){
			Q[ii1][ii2] = 0.0;
		}
		Q[ii1][max_index] = 1.0;

		/* G-S 直交化 */
		for(ii2=0; ii2<n; ii2++){
			if(zero_flags[ii2]) continue;
			tmp_product = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				Q[ii1][ii3] -= tmp_product*Q[ii2][ii3];
			}
		}

		/* G-S 直交化を再度行って精度向上 */
		for(ii2=0; ii2<n; ii2++) vect_sum[ii2] = 0.0;  /* 桁落ち防止 */
		for(ii2=0; ii2<n; ii2++){
			if(zero_flags[ii2]) continue;
			tmp_product = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				vect_sum[ii3] -= tmp_product*Q[ii2][ii3];  /* 桁落ち防止 */
			}
		}
		for(ii2=0; ii2<n; ii2++) Q[ii1][ii2] += vect_sum[ii2]; /* 桁落ち防止 */
	}

	/* 転置を元に戻す */
	transpose_overwrite_matrix(n, Q);

	RC_TRY( mm_free(zero_flags) );
	RC_TRY( mm_free(vect_sum) );

	return(NORMAL_RC);
}


RC
qr_subst (int n, double **Q, double **R, const double *b, double *x)
{
	int ii1, ii2;

	/* Q^T * b => x */
	mul_trans_matrix_vect(n, n, Q, b, x);

	/* R^-1 * b => x */
	for(ii1=n-1; ii1>=0; ii1--){
		for(ii2=ii1+1; ii2<n; ii2++){
			x[ii1] -= R[ii1][ii2] * x[ii2];
		}
		if(fabs(R[ii1][ii1]) < ABS_TOL){
			x[ii1] = 0.0;
		}else{
			x[ii1] /= R[ii1][ii1];
		}
	}

	return(NORMAL_RC);
}

