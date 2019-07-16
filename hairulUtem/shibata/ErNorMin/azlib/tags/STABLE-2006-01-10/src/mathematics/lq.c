/*********************************************************************
 * lq.c
 *
 * Copyright (C) 2004 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: lq.c 420 2005-08-12 07:02:40Z sasaoka $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"


RC
lq_solve (double **A, int n, const double *b, double *x, int active[])
{
	double **copy_A;
	double **L;
	double **inv_L;
	double **Q;
	double *fact;
	double *copy_b;
	double *scale;
	int ii1, ii2, ii3;
	double sum, max_fact;
	int max_index;

	RC_TRY( allocate2D(n, n, &copy_A) );
	RC_TRY( allocate2D(n, n, &L) );
	RC_TRY( allocate2D(n, n, &inv_L) );
	RC_TRY( allocate2D(n, n, &Q) );
	RC_NULL_CHK( fact = (double *)mm_alloc(sizeof(double)*n) );
	RC_NULL_CHK( copy_b = (double *)mm_alloc(sizeof(double)*n) );
	RC_NULL_CHK( scale = (double *)mm_alloc(sizeof(double)*n) );

	/* A, b を行毎に最大値スケーリングして copy_A, copy_b に */
	RC_TRY( copy_matrix(n, n, A, copy_A) );
	for(ii1=0; ii1<n; ii1++) copy_b[ii1] = b[ii1];
	for(ii1=0; ii1<n; ii1++){
		double scale_fact = ABS_TOL;
		for(ii2=0; ii2<n; ii2++){
			if(fabs(A[ii1][ii2]) > scale_fact){
				scale_fact = fabs(A[ii1][ii2]);
			}
		}
		for(ii2=0; ii2<n; ii2++){
			copy_A[ii1][ii2] /= scale_fact;
		}
		copy_b[ii1] /= scale_fact;
	}

	/* copy_A を列毎に最大値スケーリングしてスケーリング係数を scale に */
	for(ii1=0; ii1<n; ii1++){
		double scale_fact = ABS_TOL;
		for(ii2=0; ii2<n; ii2++){
			if(fabs(copy_A[ii2][ii1]) > scale_fact){
				scale_fact = fabs(copy_A[ii2][ii1]);
			}
		}
		for(ii2=0; ii2<n; ii2++){
			copy_A[ii2][ii1] /= scale_fact;
		}
		scale[ii1] = scale_fact;
	}

	RC_TRY( lq_factor(n, copy_A, L, Q) );

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			inv_L[ii1][ii2] = 0.0;
		}
		inv_L[ii1][ii1] = 1.0;

		for(ii2=0; ii2<n; ii2++){
			for(ii3=0; ii3<ii2; ii3++){
				inv_L[ii1][ii2] -= L[ii2][ii3] * inv_L[ii1][ii3];
			}
			if(fabs(L[ii2][ii2]) < ABS_TOL){
				inv_L[ii1][ii2] = 0.0;
			}else{
				inv_L[ii1][ii2] /= L[ii2][ii2];
			}
		}
	}
	transpose_overwrite_matrix(n, inv_L);

	if(active != NULL){
		for(ii1=0; ii1<n; ii1++) active[ii1] = 1;
	}
	for(ii1=0; ii1<n; ii1++){
		if(fabs(L[ii1][ii1]) >= ABS_TOL) continue;

		mul_vect_matrix(n, n, L[ii1], inv_L, fact);
/*		fprintf(stderr, "[%d]", ii1);
		print_vect(stderr, n, fact); */
		sum = 0.0;
		max_index = -1;
		max_fact = ABS_TOL;
		for(ii2=0; ii2<ii1; ii2++){
			sum += fact[ii2] * copy_b[ii2];
			if(fabs(fact[ii2]) > max_fact){
				max_index = ii2;
				max_fact = fabs(fact[ii2]);
			}
		}
		if(max_index >= 0){
			if(copy_b[ii1] < sum){
				copy_b[max_index] += (copy_b[ii1] - sum)/fact[max_index];
				if(active != NULL) active[max_index] = 0;
			}else{
				if(active != NULL) active[ii1] = 0;
			}
		}else{
			if(active != NULL) active[ii1] = 0;
		}
	}

	RC_TRY( lq_subst(n, L, Q, copy_b, x) );
	for(ii1=0; ii1<n; ii1++){
		x[ii1] /= scale[ii1];
	}

	RC_TRY( free2D(n, n, &copy_A) );
	RC_TRY( free2D(n, n, &L) );
	RC_TRY( free2D(n, n, &inv_L) );
	RC_TRY( free2D(n, n, &Q) );
	RC_TRY( mm_free(fact) );
	RC_TRY( mm_free(copy_b) );
	RC_TRY( mm_free(scale) );

	return(NORMAL_RC);
}


/* A[n][n] を下(左)三角L[n][n] と直交行列 Q[n][n]に分解 */
RC
lq_factor (int n, double **A, double **L, double **Q)
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

	init_matrix(n, n, L);
	RC_TRY( copy_matrix(n, n, A, Q) );

	min_norm = ABS_TOL;
	for(ii1=0; ii1<n; ii1++){
		tmp_product = sqrt( inner_product(n, Q[ii1], Q[ii1]) );
		if(tmp_product > min_norm) min_norm = tmp_product;
	}
	min_norm *= (10e+3)*REL_TOL;
	min_norm += ABS_TOL;

	for(ii1=0; ii1<n; ii1++){
		/* G-S 直交化 */
		for(ii2=0; ii2<ii1; ii2++){
			L[ii1][ii2] = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				Q[ii1][ii3] -= L[ii1][ii2]*Q[ii2][ii3];
			}
		}

		/* G-S 直交化を再度行って精度向上 */
		for(ii2=0; ii2<n; ii2++) vect_sum[ii2] = 0.0;  /* 桁落ち防止 */
		for(ii2=0; ii2<ii1; ii2++){
			tmp_product = inner_product(n, Q[ii1], Q[ii2]);
			for(ii3=0; ii3<n; ii3++){
				vect_sum[ii3] -= tmp_product*Q[ii2][ii3];  /* 桁落ち防止 */
			}
			L[ii1][ii2] += tmp_product;
		}
		for(ii2=0; ii2<n; ii2++) Q[ii1][ii2] += vect_sum[ii2]; /* 桁落ち防止 */

		L[ii1][ii1] = sqrt( inner_product(n, Q[ii1], Q[ii1]) );
		if(L[ii1][ii1] > min_norm){
			for(ii2=0; ii2<n; ii2++){
				Q[ii1][ii2] /= L[ii1][ii1];
			}
		}else{
			zero_flags[ii1] = 1;
			for(ii2=0; ii2<n; ii2++){
				Q[ii1][ii2] = 0.0;
			}
		}
	}

	/* zero_flags の立っている L の対角項をゼロにする */
	for(ii1=0; ii1<n; ii1++){
		if(zero_flags[ii1] == 0) continue;
		L[ii1][ii1] = 0.0;
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

	RC_TRY( mm_free(zero_flags) );
	RC_TRY( mm_free(vect_sum) );

	return(NORMAL_RC);
}


RC
lq_subst (int n, double **L, double **Q, const double *b, double *x)
{
	int ii1, ii2;
	double *tmp_vect;

	RC_NULL_CHK( tmp_vect = (double *)mm_alloc(sizeof(double)*n) );

	/* L^-1 * b => tmp_vect */
	for(ii1=0; ii1<n; ii1++){
		tmp_vect[ii1] = b[ii1];
		for(ii2=0; ii2<ii1; ii2++){
			tmp_vect[ii1] -= L[ii1][ii2] * tmp_vect[ii2];
		}
		if(fabs(L[ii1][ii1]) < ABS_TOL){
			tmp_vect[ii1] = 0.0;
		}else{
			tmp_vect[ii1] /= L[ii1][ii1];
		}
	}

	/* Q^T * x => x */
	mul_trans_matrix_vect(n, n, Q, tmp_vect, x);

	RC_TRY( mm_free(tmp_vect) );

	return(NORMAL_RC);
}


RC
inverse_matrix_lq (int n, double **A, double **inverseA)
{
	double **copy_A;
	double **L;
	double **Q;
	double *scale_row;
	double *scale_col;
	double *b;
	double *x;
	int ii1, ii2;

	RC_TRY( allocate2D(n, n, &copy_A) );
	RC_TRY( allocate2D(n, n, &L) );
	RC_TRY( allocate2D(n, n, &Q) );
	RC_NULL_CHK( scale_row = (double *)mm_alloc(sizeof(double)*n) );
	RC_NULL_CHK( scale_col = (double *)mm_alloc(sizeof(double)*n) );
	RC_NULL_CHK( b = (double *)mm_alloc(sizeof(double)*n) );
	RC_NULL_CHK( x = (double *)mm_alloc(sizeof(double)*n) );

	/* A を行毎に最大値スケーリングしてスケーリング係数を scale_row に */
	RC_TRY( copy_matrix(n, n, A, copy_A) );
	for(ii1=0; ii1<n; ii1++){
		double scale_fact = ABS_TOL;
		for(ii2=0; ii2<n; ii2++){
			if(fabs(A[ii1][ii2]) > scale_fact){
				scale_fact = fabs(A[ii1][ii2]);
			}
		}
		for(ii2=0; ii2<n; ii2++){
			copy_A[ii1][ii2] /= scale_fact;
		}
		scale_row[ii1] = scale_fact;
	}

	/* copy_A を列毎に最大値スケーリングしてスケーリング係数を scale_col に */
	for(ii1=0; ii1<n; ii1++){
		double scale_fact = ABS_TOL;
		for(ii2=0; ii2<n; ii2++){
			if(fabs(copy_A[ii2][ii1]) > scale_fact){
				scale_fact = fabs(copy_A[ii2][ii1]);
			}
		}
		for(ii2=0; ii2<n; ii2++){
			copy_A[ii2][ii1] /= scale_fact;
		}
		scale_col[ii1] = scale_fact;
	}

	RC_TRY( lq_factor(n, copy_A, L, Q) );

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++) b[ii2] = 0.0;
		b[ii1] = 1.0/scale_row[ii1];

		RC_TRY( lq_subst(n, L, Q, b, x) );

		for(ii2=0; ii2<n; ii2++){
			inverseA[ii1][ii2] = x[ii2]/scale_col[ii2];
		}
	}

	RC_TRY( free2D(n, n, &copy_A) );
	RC_TRY( free2D(n, n, &L) );
	RC_TRY( free2D(n, n, &Q) );
	RC_TRY( mm_free(scale_row) );
	RC_TRY( mm_free(scale_col) );
	RC_TRY( mm_free(b) );
	RC_TRY( mm_free(x) );

	return(NORMAL_RC);
}

