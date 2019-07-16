/*********************************************************************
 * math_utl.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Masanobu SUYAMA> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: math_utl.c,v 1.15 2003/12/02 10:20:47 sasaoka Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"

#define ADM_REL_ERROR (2.0*DBL_EPSILON)
#define ADM_ABS_ERROR (10.0*DBL_MIN)

static long idum1 = 1L;
static long idum2 = 123456789L;
static long iy = 0L;
static long iv[RANDOM_NTAB];

static void mul_matrix_large(int l, int m, int n,
                             double **A, double **B, double **C);
static void mul_matrix_AtBA3(int m, int n, double **A, double **B, double **C);


/* NUMERICAL RECIPES in C より */
/* 乱数列を idum で初期化する  */
void my_random_init(long idum)
{
	int j;
	long k;

	if(idum < 1){
		idum1 = 1;
	}else{
		idum1 = idum;
	}
	idum2 = 123456789;
	for(j=RANDOM_NTAB+7;j>=0;j--){
		k = idum1/RANDOM_IQ1;
		idum1 = RANDOM_IA1*(idum1 - k*RANDOM_IQ1) - k*RANDOM_IR1;
		if(idum1 < 0) idum1 += RANDOM_IM1;
		if(j < RANDOM_NTAB)iv[j] = idum1;
	}
	iy = iv[0];
}


/* 乱数列を発生させる */
double my_random(void)
{
	int j;
	long k;
	double temp;

	if(iy <= 0){
		my_random_init(1L);
	}

	k = (idum1)/RANDOM_IQ1;
	idum1 = RANDOM_IA1*(idum1 - k*RANDOM_IQ1) - k*RANDOM_IR1;
	if(idum1 < 0) idum1 += RANDOM_IM1;

	k = idum2/RANDOM_IQ2;
	idum2 = RANDOM_IA2*(idum2 - k*RANDOM_IQ2) - k*RANDOM_IR2;
	if(idum2 < 0)idum2 += RANDOM_IM2;

	j = iy/RANDOM_NDIV;
	iy = iv[j] - idum2;
	iv[j] = idum1;
	if(iy < 1)iy += RANDOM_IMM1;

	if((temp=RANDOM_AM*iy) > RANDOM_RNMX){
	 	return RANDOM_RNMX;
	}
	return temp;
}


int nearly_eq(double v1, double v2)
{
	double abs_v1_v2;
	double mean_v1_v2;
	double adm_error;

	abs_v1_v2 = fabs(v1 - v2);
	mean_v1_v2 = (fabs(v1) + fabs(v2))/2.0;
	adm_error = ADM_ABS_ERROR + (mean_v1_v2 * ADM_REL_ERROR);

	if(abs_v1_v2 < adm_error) return(1);

	return(0);
}


/* メモリの確保　〜1次元 ver.〜 */
RC allocate1D(int size, double **array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (double *)malloc(size * sizeof(double)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1] = 0.0;
	}

	return(NORMAL_RC);
}


/* メモリの確保　〜２次元 ver.〜 */
RC allocate2D(int record, int column, double ***array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (double **)malloc(record * sizeof(double *)) );

	for(ii1=0;ii1<record;ii1++){
		RC_TRY( allocate1D(column, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}


/* メモリの確保　〜３次元 ver.〜 */
RC allocate3D(int record1, int record2, int record3, double ****array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (double ***)malloc(record1 * sizeof(double **)) );

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( allocate2D(record2, record3, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}


/* メモリの確保　〜４次元 ver.〜 */
RC allocate4D(int record1, int record2,
              int record3, int record4, double *****array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (double ****)malloc(record1 * sizeof(double ***)) );

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( allocate3D(record2, record3, record4, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}

/* メモリの確保(int版)　〜1次元 ver.〜 */
RC allocate1D_i(int size, int **array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (int *)malloc(size * sizeof(int)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1] = 0;
	}

	return(NORMAL_RC);
}


/* メモリの確保(int版)　〜２次元 ver.〜 */
RC allocate2D_i(int record, int column, int ***array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (int **)malloc(record * sizeof(int *)) );

	for(ii1=0;ii1<record;ii1++){
		RC_TRY( allocate1D_i(column, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}


/* メモリの確保(int版)　〜３次元 ver.〜 */
RC allocate3D_i(int record1, int record2, int record3, int ****array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (int ***)malloc(record1 * sizeof(int **)) );

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( allocate2D_i(record2, record3, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}


/* メモリの確保(int版)　〜４次元 ver.〜 */
RC allocate4D_i(int record1, int record2,
                int record3, int record4, int *****array)
{
	int ii1;

	RC_NULL_CHK( (*array) = (int ****)malloc(record1 * sizeof(int ***)) );

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( allocate3D_i(record2, record3, record4, &((*array)[ii1])) );
	}

	return(NORMAL_RC);
}


/* メモリの解放　〜２次元 ver.〜 */
RC free2D(int record, int column, double **array)
{
	int ii1;

	for(ii1=0;ii1<record;ii1++){
		if(array[ii1] == NULL) return(ARG_ERROR_RC);
		free((void *)array[ii1]);
		array[ii1] = (double *)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}


/* メモリの解放　〜３次元 ver.〜 */
RC free3D(int record1, int record2, int record3, double ***array)
{
	int ii1;

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( free2D(record2, record3, array[ii1]) );
		array[ii1] = (double **)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}


/* メモリの解放　〜４次元 ver.〜 */
RC free4D(int record1, int record2,
          int record3, int record4, double ****array)
{
	int ii1;

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( free3D(record2, record3, record4, array[ii1]) );
		array[ii1] = (double ***)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}

/* メモリの解放(int版)　〜２次元 ver.〜 */
RC free2D_i(int record, int column, int **array)
{
	int ii1;

	for(ii1=0;ii1<record;ii1++){
		if(array[ii1] == NULL) return(ARG_ERROR_RC);
		free((void *)array[ii1]);
		array[ii1] = (int *)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}


/* メモリの解放(int版)　〜３次元 ver.〜 */
RC free3D_i(int record1, int record2, int record3, int ***array)
{
	int ii1;

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( free2D_i(record2, record3, array[ii1]) );
		array[ii1] = (int **)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}


/* メモリの解放(int版)　〜４次元 ver.〜 */
RC free4D_i(int record1, int record2,
            int record3, int record4, int ****array)
{
	int ii1;

	for(ii1=0;ii1<record1;ii1++){
		RC_TRY( free3D_i(record2, record3, record4, array[ii1]) );
		array[ii1] = (int ***)NULL;
	}

	free((void *)array);

	return(NORMAL_RC);
}

static void mul_matrix_large(int l, int m, int n,
                             double **A, double **B, double **C)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<l; ii1++){
		double *C_ii1 = C[ii1];

		for(ii2=0; ii2<n; ii2++){
			C_ii1[ii2] = 0.0;
		}
	}

	for(ii1=0; ii1<l-1; ii1+=2){
		double *C_ii1   = C[ii1];
		double *C_ii1p1 = C[ii1+1];

		for(ii2=0; ii2<m-1; ii2+=2){
			double A_ii1_ii2   = A[ii1][ii2];
			double A_ii1_ii2p1 = A[ii1][ii2+1];
			double A_ii1p1_ii2   = A[ii1+1][ii2];
			double A_ii1p1_ii2p1 = A[ii1+1][ii2+1];
			double *B_ii2 = B[ii2];
			double *B_ii2p1 = B[ii2+1];

			for(ii3=0; ii3<n; ii3++){
				double B_ii2_ii3   = B_ii2[ii3];
				double B_ii2p1_ii3 = B_ii2p1[ii3];

				C_ii1[ii3]   += A_ii1_ii2    *B_ii2_ii3
				              + A_ii1_ii2p1  *B_ii2p1_ii3;
				C_ii1p1[ii3] += A_ii1p1_ii2  *B_ii2_ii3
				              + A_ii1p1_ii2p1*B_ii2p1_ii3;
			}
		}
		for(; ii2<m; ii2++){
			double A_ii1_ii2 = A[ii1][ii2];
			double A_ii1p1_ii2   = A[ii1+1][ii2];
			double *B_ii2 = B[ii2];

			for(ii3=0; ii3<n; ii3++){
				double B_ii2_ii3   = B_ii2[ii3];

				C_ii1[ii3] += A_ii1_ii2 * B_ii2_ii3;
				C_ii1p1[ii3] += A_ii1p1_ii2 * B_ii2_ii3;
			}
		}
	}
	for(; ii1<l; ii1++){
		double *C_ii1 = C[ii1];
		double A_ii1_0 = A[ii1][0];
		double *B_0 = B[0];

		for(ii3=0; ii3<n; ii3++){
			C_ii1[ii3] = A_ii1_0 * B_0[ii3];
		}

		for(ii2=1; ii2<m; ii2++){
			double A_ii1_ii2 = A[ii1][ii2];
			double *B_ii2 = B[ii2];

			for(ii3=0; ii3<n; ii3++){
				C_ii1[ii3] += A_ii1_ii2 * B_ii2[ii3];
			}
		}
	}
}


/* multiplication of matrices A[l][m] * B[m][n] => C[l][n] */
/* 行列の掛け算                                            */
void mul_matrix(int l, int m, int n, double **A, double **B, double **C)
{

	if((l <= 0)||(m <= 0)||(n <= 0)) return;

	if(l*m*n >= 100){
		mul_matrix_large(l, m, n, A, B, C);
	}else{
		int ii1, ii2, ii3;

		for(ii1=0; ii1<l; ii1++){
			double *C_ii1 = C[ii1];
			double A_ii1_0 = A[ii1][0];
			double *B_0 = B[0];

			for(ii3=0; ii3<n; ii3++){
				C_ii1[ii3] = A_ii1_0 * B_0[ii3];
			}

			for(ii2=1; ii2<m; ii2++){
				double A_ii1_ii2 = A[ii1][ii2];
				double *B_ii2 = B[ii2];

				for(ii3=0; ii3<n; ii3++){
					C_ii1[ii3] += A_ii1_ii2 * B_ii2[ii3];
				}
			}
		}
	}
}


/* A[m][n] * X[n] => Y[m] */
void mul_matrix_vect(int m, int n, double **A, const double *X, double *Y)
{
	int ii1, ii2;
	long double sum;

	for(ii1=0;ii1<m;ii1++){
		sum = 0.0;

		for(ii2=0;ii2<n;ii2++){
			sum += A[ii1][ii2] * X[ii2];
		}
		Y[ii1] = (double)sum;
	}
}

/* A^T[n][m] * X[m] => Y[n] */
void mul_trans_matrix_vect(int m, int n, double **A, const double *X, double *Y)
{
	int ii1, ii2;
	long double sum;

	for(ii1=0;ii1<n;ii1++){
		sum = 0.0;

		for(ii2=0;ii2<m;ii2++){
			sum += A[ii2][ii1] * X[ii2];
		}
		Y[ii1] = (double)sum;
	}
}


/* X^T[m] * A[m][n] => Y[n] */
void mul_vect_matrix(int m, int n, const double *X, double **A, double *Y)
{
	int ii1, ii2;
	long double sum;

	for(ii1=0;ii1<n;ii1++){
		sum = 0.0;

		for(ii2=0;ii2<m;ii2++){
			sum += X[ii2] * A[ii2][ii1];
		}
		Y[ii1] = (double)sum;
	}
}


/* 行列の行列式 */
double determinant(int dim, double **matrix)
{
	if(dim == 1){
		return(matrix[0][0]);
	}else if(dim == 2){
		return(determinant2(matrix));
	}else if(dim == 3){
		return(determinant3(matrix));
	}else{
		return(0.0);
	}
}


/* 2 x 2 行列の行列式 */
double determinant2(double **matrix)
{
	return(matrix[0][0] * matrix[1][1]
	     - matrix[0][1] * matrix[1][0]);
}


/* 3 x 3 行列の行列式 */
double determinant3(double **matrix)
{
	return(matrix[0][0] * matrix[1][1] * matrix[2][2]
	     + matrix[0][1] * matrix[1][2] * matrix[2][0]
		 + matrix[0][2] * matrix[1][0] * matrix[2][1]
		 - matrix[0][2] * matrix[1][1] * matrix[2][0]
		 - matrix[0][1] * matrix[1][0] * matrix[2][2]
		 - matrix[0][0] * matrix[1][2] * matrix[2][1]);
}


/* A の逆マトリックスを inverseA に代入    */
/* dim=2 ...2x2, dim=3 ...3x3 dim=n ...nxn */
RC inverse_matrix(int dim, double **A, double **inverseA)
{
	if(dim < 1) return(ARG_ERROR_RC);

	if(dim == 1){
		if(nearly_eq(A[0][0], 0.0)) return(CAL_ERROR_RC);
		inverseA[0][0] = 1.0/A[0][0];
	}else if(dim == 2){
		double det;

		det = determinant2(A);
		if(nearly_eq(det, 0.0)) return(CAL_ERROR_RC);
		inverseA[0][0] = A[1][1] / det;
		inverseA[0][1] = - A[0][1] / det;
		inverseA[1][0] = - A[1][0] / det;
		inverseA[1][1] = A[0][0] / det;
	}else if(dim == 3){
		double det;

		det = determinant3(A);
		if(nearly_eq(det, 0.0)) return(CAL_ERROR_RC);
		inverseA[0][0] = (- A[1][2] * A[2][1] + A[1][1] * A[2][2])/det;
		inverseA[0][1] = (  A[0][2] * A[2][1] - A[0][1] * A[2][2])/det;
		inverseA[0][2] = (- A[0][2] * A[1][1] + A[0][1] * A[1][2])/det;
		inverseA[1][0] = (  A[1][2] * A[2][0] - A[1][0] * A[2][2])/det;
		inverseA[1][1] = (- A[0][2] * A[2][0] + A[0][0] * A[2][2])/det;
		inverseA[1][2] = (  A[0][2] * A[1][0] - A[0][0] * A[1][2])/det;
		inverseA[2][0] = (- A[1][1] * A[2][0] + A[1][0] * A[2][1])/det;
		inverseA[2][1] = (  A[0][1] * A[2][0] - A[0][0] * A[2][1])/det;
		inverseA[2][2] = (- A[0][1] * A[1][0] + A[0][0] * A[1][1])/det;
	}else{
		int ii1, ii2;
		double dum;
		double *col;
		double **LU;
		int *index;

		col = (double *)malloc(dim * sizeof(double));
		index = (int *)malloc(dim * sizeof(int));
		if((col == NULL)||(index == NULL)){
			fprintf(stderr, "malloc() < ");
			return(ALLOC_ERROR_RC);
		}
		RC_TRY( allocate2D(dim, dim, &LU) );
		RC_TRY( copy_matrix(dim, dim, A, LU) );

		RC_TRY( lu_decomp(LU, dim, index, &dum) );

		for(ii1=0; ii1<dim; ii1++){
			for(ii2=0; ii2<dim; ii2++){
				col[ii2] = 0.0;
			}
			col[ii1] = 1.0;

			RC_TRY( lu_subst(LU, dim, index, col) );
			for(ii2=0; ii2<dim; ii2++){
				inverseA[ii2][ii1] = col[ii2];
			}
		}

		free(col);
		free(index);
		RC_TRY( free2D(dim, dim, LU) );
	}

	return(NORMAL_RC);
}


/* wa * A[m][n] + wb * B[m][n] => C[m][n] */
/* A == C or B == C でも構わない          */
RC add_matrix(int m, int n, double wa, double **A, double wb, double **B,
              double **C)
{
	int ii1 ,ii2;

	if((m <= 0)||(n <= 0)||(A == NULL)||(B == NULL)||(C == NULL)){
		return(ARG_ERROR_RC);
	}
	
	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n ; ii2++){
			C[ii1][ii2] = (wa*(A[ii1][ii2])) + (wb*(B[ii1][ii2]));
		}
	}

	return(NORMAL_RC);
}


/* A[m][n] を B[m][n] にコピー */
RC copy_matrix(int m, int n, double **A, double **B)
{
	int ii1 ,ii2;

	if((m <= 0)||(n <= 0)||(A == NULL)||(B == NULL)){
		return(ARG_ERROR_RC);
	}
	
	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n ; ii2++){
			B[ii1][ii2] = A[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


/* Ax = b を解く                     */
/* x と b は同じポインタでも構わない */
/* A の内容は変化しない              */
RC lu_solve(double **A, int n, const double *b, double *x)
{
	double **LU;
	int *index;
	double d;
	int ii1;

	RC_TRY( allocate2D(n, n, &LU) );
	RC_NULL_CHK( index = malloc(n*sizeof(int)) );

	RC_TRY( copy_matrix(n, n, A, LU) );
	RC_TRY( lu_decomp(LU, n, index, &d) );
	for(ii1=0; ii1<n; ii1++){
		x[ii1] = b[ii1];
	}
	RC_TRY( lu_subst(LU, n, index, x) );

	RC_TRY( free2D(n, n, LU) );
	free(index);

	return(NORMAL_RC);
}


/* 行列 A[n][n] を LU 分解して A に上書き        */
/* index[n] : ピボット選択による行交換の情報     */
/* d : 行交換の回数が偶数なら 1.0, 奇数なら -1.0 */
RC lu_decomp(double **A, int n, int *index, double *d)
{
	int ii1, ii2, ii3;
	double *scale_factor;
	double max;
	double tmp;
	int max_index;

	if((n <= 0)||(A == NULL)||(index == NULL)||(d == NULL)){
		return(ARG_ERROR_RC);
	}

	scale_factor = (double *)malloc(n * sizeof(double));
	if(scale_factor == NULL){
		fprintf(stderr, "malloc() < ");
		return(ALLOC_ERROR_RC);
	}

	/* scale_factor を設定 */
	for(ii1=0; ii1<n; ii1++){
		max = 0.0;
		for(ii2=0; ii2<n; ii2++){
			tmp = fabs(A[ii1][ii2]);
			if(tmp > max) max = tmp;
		}
		if( nearly_eq(max, 0.0) ) return(CAL_ERROR_RC);
		scale_factor[ii1] = 1.0/max;
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
		max_index = -1;
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

		if(max_index < 0) return(CAL_ERROR_RC);
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

	free(scale_factor);

	return(NORMAL_RC);
}


/*
 * 前進代入、後退代入を行なう
 * A[n][n], index[n] : lu_decomp で得た値
 * B[n]              : 右辺ベクトル -> 解ベクトル
 * B[] に多くのゼロが含まれる場合(逆マトリックスの計算時)を想定した高速化
 */
RC lu_subst(double **A, int n, int *index, double *B)
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
		if((index[ii1] < 0)||(index[ii1] >= n)){
			fprintf(stderr, "index[ii1] < ");
			return(ARG_ERROR_RC);
		}
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
		if(nearly_eq(A[ii1][ii1], 0.0)) return(CAL_ERROR_RC);
		B[ii1] /= A[ii1][ii1];
	}

	return(NORMAL_RC);
}


/* 行列 A (row * col) を転置して B に代入 */
void transpose_matrix(int row, int col, double **A, double **B)
{
	int ii1, ii2;

	for(ii1=0; ii1<col; ii1++){
		for(ii2=0; ii2<row; ii2++){
			B[ii1][ii2] = A[ii2][ii1];
		}
	}
}


/* 正方行列 A[n][n] を転置して上書き */
void transpose_overwrite_matrix(int n, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<n; ii1++){
		for(ii2=ii1+1; ii2<n; ii2++){
			double tmp = A[ii1][ii2];
			A[ii1][ii2] = A[ii2][ii1];
			A[ii2][ii1] = tmp;
		}
	}
}


/* 行列 A (row * col) を 0.0 で初期化 */
void init_matrix(int row, int col, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<row; ii1++){
		for(ii2=0; ii2<col; ii2++){
			A[ii1][ii2] = 0.0;
		}
	}
}


/* 行列 A (n * n) を単位行列として初期化 */
void unit_matrix(int n, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			A[ii1][ii2] = 0.0;
		}
		A[ii1][ii1] = 1.0;
	}
}


TRANSLATION2D init_translation2d(void)
{
	TRANSLATION2D ret = {0.0, 0.0};
	return(ret);
}


TRANSLATION3D init_translation3d(void)
{
	TRANSLATION3D ret = {0.0, 0.0, 0.0};
	return(ret);
}


ROTATION3D init_rotation3d(void)
{
	ROTATION3D ret = {0.0, 0.0, 0.0};
	return(ret);
}


TRANS_ROTATION3D init_trans_rotation3d(void)
{
	TRANS_ROTATION3D ret = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	return(ret);
}


I_TRANS_ROTATION3D init_i_trans_rotation3d(void) 
{
	I_TRANS_ROTATION3D ret = {0, 0, 0, 0, 0, 0};
	return(ret);
}


/* p1, p2 の中央に位置する点 〜TRANSLATION2D ver.〜 */
TRANSLATION2D mid_point2d(TRANSLATION2D p1, TRANSLATION2D p2)
{
	TRANSLATION2D ret;

	ret.x = (p1.x + p2.x)/2.0;
	ret.y = (p1.y + p2.y)/2.0;

	return(ret);
}


/* p1, p2 の中央に位置する点 〜TRANSLATION3D ver.〜 */
TRANSLATION3D mid_point3d(TRANSLATION3D p1, TRANSLATION3D p2)
{
	TRANSLATION3D ret;

	ret.x = (p1.x + p2.x)/2.0;
	ret.y = (p1.y + p2.y)/2.0;
	ret.z = (p1.z + p2.z)/2.0;

	return(ret);
}


/* ベクトルの絶対値 〜TRANSLATION2D ver.〜 */
double abs_translation2d(TRANSLATION2D v)
{
	return(sqrt(v.x*v.x + v.y*v.y));
}


/* ベクトルの絶対値 〜TRANSLATION3D ver.〜 */
double abs_translation3d(TRANSLATION3D v)
{
	return(sqrt(v.x*v.x + v.y*v.y + v.z*v.z));
}


/* ROTATION3Dの絶対値 */
double abs_rotation3d(ROTATION3D v)
{
	return(sqrt(v.yz*v.yz + v.zx*v.zx + v.xy*v.xy));
}


/* p1, p2 間の距離 〜TRANSLATION2D ver.〜 */
double dist_point2d(TRANSLATION2D p1, TRANSLATION2D p2)
{
	double ret;

	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	ret = sqrt(dx*dx + dy*dy);

	return(ret);
}


/* p1, p2 間の距離 〜TRANSLATION3D ver.〜 */
double dist_point3d(TRANSLATION3D p1, TRANSLATION3D p2)
{
	double ret;

	double dx = p1.x - p2.x;
	double dy = p1.y - p2.y;
	double dz = p1.z - p2.z;
	ret = sqrt(dx*dx + dy*dy + dz*dz);

	return(ret);
}


/* ベクトルの引き算 (v1 - v2) 〜TRANSLATION2D ver.〜 */
TRANSLATION2D sub_translation2d(TRANSLATION2D v1, TRANSLATION2D v2)
{
	TRANSLATION2D v;

	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;

	return(v);
}


/* ベクトルの引き算 (v1 - v2) 〜TRANSLATION3D ver.〜 */
TRANSLATION3D sub_translation3d(TRANSLATION3D v1, TRANSLATION3D v2)
{
	TRANSLATION3D v;

	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	v.z = v1.z - v2.z;

	return(v);
}


/* ベクトルの足し算 (v1 + v2) 〜TRANSLATION2D ver.〜 */
TRANSLATION2D add_translation2d(TRANSLATION2D v1, TRANSLATION2D v2)
{
	TRANSLATION2D v;

	v.x = v2.x + v1.x;
	v.y = v2.y + v1.y;

	return(v);
}


/* ベクトルの足し算 (v1 + v2) 〜TRANSLATION3D ver.〜 */
TRANSLATION3D add_translation3d(TRANSLATION3D v1, TRANSLATION3D v2)
{
	TRANSLATION3D v;

	v.x = v2.x + v1.x;
	v.y = v2.y + v1.y;
	v.z = v2.z + v1.z;

	return(v);
}


TRANS_ROTATION3D add_trans_rotation3d(TRANS_ROTATION3D v1, TRANS_ROTATION3D v2)
{
	TRANS_ROTATION3D v;

	v.x = v2.x + v1.x;
	v.y = v2.y + v1.y;
	v.z = v2.z + v1.z;
	v.yz = v2.yz + v1.yz;
	v.zx = v2.zx + v1.zx;
	v.xy = v2.xy + v1.xy;

	return(v);
}


TRANS_ROTATION3D ave_trans_rotation3d(TRANS_ROTATION3D v1, TRANS_ROTATION3D v2)
{
	TRANS_ROTATION3D v;

	v.x = (v2.x + v1.x)/2.0;
	v.y = (v2.y + v1.y)/2.0;
	v.z = (v2.z + v1.z)/2.0;
	v.yz = (v2.yz + v1.yz)/2.0;
	v.zx = (v2.zx + v1.zx)/2.0;
	v.xy = (v2.xy + v1.xy)/2.0;

	return(v);
}


/* スカラー三重積 v1 (dot) v2 (cross) v3 */
double scalar_triple_product3d(TRANSLATION3D v1,
		                       TRANSLATION3D v2, TRANSLATION3D v3)
{
	double s;
	TRANSLATION3D v;

	v = outer_product3d(v2, v3);
	s = inner_product3d(v1, v);

	return(s);
}


/* ベクトルの内積 v1 (dot) v2 〜TRANSLATION2D ver.〜 */
double inner_product2d(TRANSLATION2D v1, TRANSLATION2D v2)
{
	return(v1.x * v2.x + v1.y * v2.y);
}


/* ベクトルの内積 v1 (dot) v2 〜TRANSLATION3D ver.〜 */
double inner_product3d(TRANSLATION3D v1, TRANSLATION3D v2)
{
	return(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z);
}


/* ベクトルの外積 v1 (cross) v2 〜TRANSLATION2D ver.〜 */
double outer_product2d(TRANSLATION2D v1, TRANSLATION2D v2)
{
	return((v1.x * v2.y) - (v1.y * v2.x));
}


/* ベクトルの外積 v1 (cross) v2 〜TRANSLATION3D ver.〜 */
TRANSLATION3D outer_product3d(TRANSLATION3D v1, TRANSLATION3D v2)
{
	TRANSLATION3D v;

	v.x = (v1.y * v2.z) - (v1.z * v2.y);
	v.y = (v1.z * v2.x) - (v1.x * v2.z);
	v.z = (v1.x * v2.y) - (v1.y * v2.x);

	return(v);
}


/* ベクトルのテンソル積 v1 * v2 => ans[][] */
double **tensor_product3d(TRANSLATION3D v1, TRANSLATION3D v2, double **ans)
{
	ans[0][0] = v1.x * v2.x;
	ans[0][1] = v1.x * v2.y;
	ans[0][2] = v1.x * v2.z;
	ans[1][0] = v1.y * v2.x;
	ans[1][1] = v1.y * v2.y;
	ans[1][2] = v1.y * v2.z;
	ans[2][0] = v1.z * v2.x;
	ans[2][1] = v1.z * v2.y;
	ans[2][2] = v1.z * v2.z;

	return(ans);
}


RC allocate_operation_matrix(double ***matrix)
{
	RC_TRY( allocate2D(4, 4, matrix) );
	RC_TRY( init_operation_matrix(*matrix) );

	return(NORMAL_RC);
}


void free_operation_matrix(double **matrix){
	free2D(4, 4, matrix);
}


RC init_operation_matrix(double **matrix)
{
	int ii1, ii2;

	for(ii1=0; ii1<4; ii1++){
		for(ii2=0; ii2<4; ii2++){
			matrix[ii1][ii2] = 0.0;
		}
		matrix[ii1][ii1] = 1.0;
	}

	return(NORMAL_RC);
}


RC set_operation_matrix_pxpy(double fact, double **mat)
{
	mat[0][0] = fact; mat[0][1] = 0.0;  mat[0][2] = 0.0;  mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = fact; mat[1][2] = 0.0;  mat[1][3] = 0.0;
	mat[2][0] = 0.0;  mat[2][1] = 0.0;  mat[2][2] = fact; mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC set_operation_matrix_mzpy(double fact, double **mat)
{
	mat[0][0] = 0.0;  mat[0][1] = 0.0;  mat[0][2] =-fact; mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = fact; mat[1][2] = 0.0;  mat[1][3] = 0.0;
	mat[2][0] = fact; mat[2][1] = 0.0;  mat[2][2] = 0.0;  mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC set_operation_matrix_pxmz(double fact, double **mat)
{
	mat[0][0] = fact; mat[0][1] = 0.0;  mat[0][2] = 0.0;  mat[0][3] = 0.0;
	mat[1][0] = 0.0;  mat[1][1] = 0.0;  mat[1][2] =-fact; mat[1][3] = 0.0;
	mat[2][0] = 0.0;  mat[2][1] = fact; mat[2][2] = 0.0;  mat[2][3] = 0.0;
	mat[3][0] = 0.0;  mat[3][1] = 0.0;  mat[3][2] = 0.0;  mat[3][3] = 1.0;

	return(NORMAL_RC);
}


RC shift_operation_matrix(double dx, double dy, double dz, double **mat)
{
	mat[0][0] += dx*mat[3][0];
	mat[0][1] += dx*mat[3][1];
	mat[0][2] += dx*mat[3][2];
	mat[0][3] += dx*mat[3][3];
	mat[1][0] += dy*mat[3][0];
	mat[1][1] += dy*mat[3][1];
	mat[1][2] += dy*mat[3][2];
	mat[1][3] += dy*mat[3][3];
	mat[2][0] += dz*mat[3][0];
	mat[2][1] += dz*mat[3][1];
	mat[2][2] += dz*mat[3][2];
	mat[2][3] += dz*mat[3][3];

	return(NORMAL_RC);
}

RC scale_operation_matrix(double x, double y, double z, double **mat)
{
	int ii1;
		 
	for(ii1=0;ii1<4;ii1++){
		mat[0][ii1] *= x; 
		mat[1][ii1] *= y;
		mat[2][ii1] *= z;
	}

	return(NORMAL_RC);
}

RC rotation_operation_matrix_x(double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[1][ii1];
		m[1][ii1] = mat[2][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[1][0] = c*m[0][0] - s*m[1][0];
	mat[1][1] = c*m[0][1] - s*m[1][1];
	mat[1][2] = c*m[0][2] - s*m[1][2];
	mat[1][3] = c*m[0][3] - s*m[1][3];
	mat[2][0] = s*m[0][0] + c*m[1][0];
	mat[2][1] = s*m[0][1] + c*m[1][1];
	mat[2][2] = s*m[0][2] + c*m[1][2];
	mat[2][3] = s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}

RC rotation_operation_matrix_y(double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[0][ii1];
		m[1][ii1] = mat[2][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[0][0] = c*m[0][0] + s*m[1][0];
	mat[0][1] = c*m[0][1] + s*m[1][1];
	mat[0][2] = c*m[0][2] + s*m[1][2];
	mat[0][3] = c*m[0][3] + s*m[1][3];
	mat[2][0] = -s*m[0][0] + c*m[1][0];
	mat[2][1] = -s*m[0][1] + c*m[1][1];
	mat[2][2] = -s*m[0][2] + c*m[1][2];
	mat[2][3] = -s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}

RC rotation_operation_matrix_z(double rad, double **mat)
{
	double m[2][4];
	double c,s;
	int ii1;

	for(ii1=0;ii1<4;ii1++){
		m[0][ii1] = mat[0][ii1];
		m[1][ii1] = mat[1][ii1];
	}
	c = cos(rad);
	s = sin(rad);

	mat[0][0] = c*m[0][0] - s*m[1][0];
	mat[0][1] = c*m[0][1] - s*m[1][1];
	mat[0][2] = c*m[0][2] - s*m[1][2];
	mat[0][3] = c*m[0][3] - s*m[1][3];
	mat[1][0] = s*m[0][0] + c*m[1][0];
	mat[1][1] = s*m[0][1] + c*m[1][1];
	mat[1][2] = s*m[0][2] + c*m[1][2];
	mat[1][3] = s*m[0][3] + c*m[1][3];

	return(NORMAL_RC);
}


/* 変換マトリックス B[4][4] に A[4][4] を作用させる */
RC mul_operation_matrix(double **A, double **B)
{
	double m[4][4];
	int ii1, ii2;

	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<4;ii2++){
			m[ii1][ii2] = B[ii1][ii2];
		}
	}

	for(ii1=0;ii1<4;ii1++){
		for(ii2=0;ii2<4;ii2++){
			B[ii1][ii2] = A[ii1][0] * m[0][ii2] + A[ii1][1] * m[1][ii2]
			            + A[ii1][2] * m[2][ii2] + A[ii1][3] * m[3][ii2];
		}
	}

	return(NORMAL_RC);
}


/* ベクトル v に変換マトリックス matrix[4][4] を作用させる */
TRANSLATION3D operation3d(double **matrix, TRANSLATION3D v)
{
	TRANSLATION3D ret;

	ret.x = matrix[0][0]*v.x + matrix[0][1]*v.y
	      + matrix[0][2]*v.z + matrix[0][3];
	ret.y = matrix[1][0]*v.x + matrix[1][1]*v.y
	      + matrix[1][2]*v.z + matrix[1][3];
	ret.z = matrix[2][0]*v.x + matrix[2][1]*v.y
	      + matrix[2][2]*v.z + matrix[2][3];

	return(ret);
}


/* TRANSLATION2D の出力 */
RC print_translation2d(FILE *fp, TRANSLATION2D v)
{
	fprintf(fp, "  x: %10.3e, y: %10.3e\n", v.x, v.y);

	return(NORMAL_RC);
}


/* TRANSLATION3D の出力 */
RC print_translation3d(FILE *fp, TRANSLATION3D v)
{
	fprintf(fp, "  x: %10.3e, y: %10.3e, z: %10.3e\n", v.x, v.y, v.z);

	return(NORMAL_RC);
}


/* ROTATION3D の出力 */
RC print_rotation3d(FILE *fp, ROTATION3D v)
{
	fprintf(fp, "  yz:%10.3e, zx:%10.3e, xy:%10.3e\n", v.yz, v.zx, v.xy);

	return(NORMAL_RC);
}


/* TRANS_ROTATION3D の出力 */
RC print_trans_rotation3d(FILE *fp, TRANS_ROTATION3D v)
{
	fprintf(fp, "  x: %10.3e, y: %10.3e, z: %10.3e\n"
	            "  yz:%10.3e, zx:%10.3e, xy:%10.3e\n",
	             v.x, v.y, v.z, v.yz, v.zx, v.xy);

	return(NORMAL_RC);
}


/* I_TRANS_ROTATION3D の出力 */
RC print_i_trans_rotation3d(FILE *fp, I_TRANS_ROTATION3D v)
{
	fprintf(fp, "  ix: %9d, iy: %9d, iz: %9d\n"
	            "  iyz:%9d, izx:%9d, ixy:%9d\n",
	             v.ix, v.iy, v.iz, v.iyz, v.izx, v.ixy);

	return(NORMAL_RC);
}


/* v1[m] と v2[n] のテンソル積 => A[m][n] */
RC tensor_product(int m, int n, const double *v1,
                                const double *v2, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n; ii2++){
			A[ii1][ii2] = v1[ii1] * v2[ii2];
		}
	}

	return(NORMAL_RC);
}


/* v1[n] と v2[n] の内積 */
double inner_product(int n, const double *v1, const double *v2)
{
	long double ret = 0.0;
	int ii1;

	for(ii1=0; ii1<n; ii1++){
		ret += v1[ii1]*v2[ii1];
	}

	return((double)ret);
}


/* A[m][n]^T * B[m][m] * A[m][n] = C[n][n] */
void mul_matrix_AtBA(int m, int n, double **A, double **B, double **C)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<n; ii1++){
		double *C_ii1 = C[ii1];

		for(ii2=0; ii2<n; ii2++){
			C_ii1[ii2] = 0.0;
		}
	}

	if(m%3 == 0){
		mul_matrix_AtBA3(m, n, A, B, C);
		return;
	}

	for(ii1=0; ii1<n; ii1++){
		double *C_ii1 = C[ii1];

		for(ii2=0; ii2<m; ii2++){
			double sum = 0.0;
			double *A_ii2 = A[ii2];

			for(ii3=0; ii3<m; ii3++){
				sum += A[ii3][ii1] * B[ii3][ii2];
			}
			for(ii3=0; ii3<n; ii3++){
				C_ii1[ii3] += sum * A_ii2[ii3];
			}
		}
	}
}


/* m%3 == 0 */
static void mul_matrix_AtBA3(int m, int n, double **A, double **B, double **C)
{
	int ii1, ii2, ii3;

	for(ii1=0; ii1<n; ii1++){
		double *C_ii1 = C[ii1];

		for(ii2=0; ii2<m; ii2+=3){
			double sum0 = 0.0;
			double sum1 = 0.0;
			double sum2 = 0.0;
			double *A_ii2 = A[ii2];
			double *A_ii2p1 = A[ii2+1];
			double *A_ii2p2 = A[ii2+2];

			for(ii3=0; ii3<m; ii3+=3){
				double A_ii3_ii1   = A[ii3][ii1];
				double A_ii3p1_ii1 = A[ii3+1][ii1];
				double A_ii3p2_ii1 = A[ii3+2][ii1];
				double *B_ii3 = B[ii3];
				double *B_ii3p1 = B[ii3+1];
				double *B_ii3p2 = B[ii3+2];

				sum0 += A_ii3_ii1  *B_ii3[ii2]
				      + A_ii3p1_ii1*B_ii3p1[ii2]
				      + A_ii3p2_ii1*B_ii3p2[ii2];
				sum1 += A_ii3_ii1  *B_ii3[ii2+1]
				      + A_ii3p1_ii1*B_ii3p1[ii2+1]
				      + A_ii3p2_ii1*B_ii3p2[ii2+1];
				sum2 += A_ii3_ii1  *B_ii3[ii2+2]
				      + A_ii3p1_ii1*B_ii3p1[ii2+2]
				      + A_ii3p2_ii1*B_ii3p2[ii2+2];
			}
			for(ii3=0; ii3<n; ii3++){
				C_ii1[ii3] += sum0*A_ii2[ii3]
				            + sum1*A_ii2p1[ii3]
				            + sum2*A_ii2p2[ii3];
			}
		}
	}
}


/* print A[m][n] */
void print_matrix(FILE *fp, int m, int n, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n; ii2++){
			fprintf(fp, " %10.3e", A[ii1][ii2]);
		}
		fprintf(fp, "\n");
	}
}


/* fa * A 〜TRANSLATION2D ver.〜 */
TRANSLATION2D mul_scalar_translation2d(double fa, TRANSLATION2D A)
{
	TRANSLATION2D ret;

	ret.x = fa * A.x;
	ret.y = fa * A.y;

	return(ret);
}


/* fa * A 〜TRANSLATION3D ver.〜 */
TRANSLATION3D mul_scalar_translation3d(double fa, TRANSLATION3D A)
{
	TRANSLATION3D ret;

	ret.x = fa * A.x;
	ret.y = fa * A.y;
	ret.z = fa * A.z;

	return(ret);
}


/* fa * A */
TRANS_ROTATION3D mul_scalar_trans_rotation3d(double fa, TRANS_ROTATION3D A)
{
	TRANS_ROTATION3D ret;

	ret.x = fa * A.x;
	ret.y = fa * A.y;
	ret.z = fa * A.z;
	ret.yz = fa * A.yz;
	ret.zx = fa * A.zx;
	ret.xy = fa * A.xy;

	return(ret);
}


TRANSLATION3D transrot2trans(TRANS_ROTATION3D A)
{
	TRANSLATION3D ret;

	ret.x = A.x;
	ret.y = A.y;
	ret.z = A.z;

	return(ret);
}


ROTATION3D transrot2rot(TRANS_ROTATION3D A)
{
	ROTATION3D ret;

	ret.yz = A.yz;
	ret.zx = A.zx;
	ret.xy = A.xy;

	return(ret);
}


/* 配列を TRANSLATION3D に変換 */
TRANSLATION3D array2translation3d(double *array)
{
	TRANSLATION3D ret;

	ret.x = array[0];
	ret.y = array[1];
	ret.z = array[2];

	return(ret);
}


RC principal3d(TRANS_ROTATION3D ss_v, double pr_ss[],
               double pr_dir0[], double pr_dir1[], double pr_dir2[])
{
	static int alloc_flag = 0;
	static double **v_matrix = NULL;
	static double **dir_matrix = NULL;
	int ii1;

	if(alloc_flag == 0){
		RC_TRY( allocate2D(3, 3, &v_matrix) );
		RC_TRY( allocate2D(3, 3, &dir_matrix) );
		alloc_flag = 1;
	}

	v_matrix[0][0] = ss_v.x;
	v_matrix[1][1] = ss_v.y;
	v_matrix[2][2] = ss_v.z;
	v_matrix[0][1] = v_matrix[1][0] = ss_v.xy;
	v_matrix[1][2] = v_matrix[2][1] = ss_v.yz;
	v_matrix[2][0] = v_matrix[0][2] = ss_v.zx;

	/* 固有値 -> pr_ss */
	/* 固有ベクトル -> dir_matrix の各列 */
	RC_TRY( eigen_jacobi_trans(3, v_matrix, pr_ss, dir_matrix) );

	for(ii1=0; ii1<3; ii1++){
		if(pr_dir0 != NULL) pr_dir0[ii1] = dir_matrix[ii1][0];
		if(pr_dir1 != NULL) pr_dir1[ii1] = dir_matrix[ii1][1];
		if(pr_dir2 != NULL) pr_dir2[ii1] = dir_matrix[ii1][2];
	}

	return(NORMAL_RC);

}


#define J_TRANS_MAX 100

/* Jacobi 変換による対称行列a[n][n]の固有値e_vals[n]、*/
/* 固有ベクトルe_vects[n][n]の計算 */
RC eigen_jacobi_trans(int n, double **a, double *e_vals, double **e_vects)
{
	int ii1, ii2, ii3, trans_count;
	double sum, sum_d, thresh;
	double *add_d;      /* 対角要素に加算する数値 */

	if(n < 1){
		return(ARG_ERROR_RC);
	}

	unit_matrix(n, e_vects);

	add_d = (double *)malloc(n * sizeof(double));
	if(add_d == NULL){
		fprintf(stderr, "malloc < ");
		return(ALLOC_ERROR_RC);
	}
	for(ii1=0; ii1<n; ii1++){
		e_vals[ii1] = a[ii1][ii1];
		add_d[ii1] = 0.0;
	}

	trans_count = 0;
	while(1){
		if(trans_count > J_TRANS_MAX){
			fprintf(stderr, "trans_count < ");
			return(CAL_ERROR_RC);
		}
		trans_count++;

		sum = 0.0;
		sum_d = 0.0;
		for(ii1=0; ii1<n; ii1++){
			for(ii2=ii1+1; ii2<n; ii2++){
				sum += fabs(a[ii1][ii2]);
			}
			sum_d += fabs(a[ii1][ii1]);
		}

		if(sum <= (10.0*DBL_MIN)) break;   /* 絶対収束 */
		
		if((sum/(double)n) <= ((10.0*DBL_EPSILON)*sum_d)){   /* 収束判定 */
			break;
		}

		if(trans_count <= 3){
			thresh = ((0.20/3.0)*(4-trans_count))*sum/((double)(n * n));
		}else{
			thresh = 0.0;
		}

		for(ii1=0; ii1<n; ii1++){
			for(ii2=ii1+1; ii2<n; ii2++){
				double g = 100.0*fabs(a[ii1][ii2]);
				double h, t, c, s, tau, t_a;

				if(fabs(a[ii1][ii2]) <= thresh){
					continue;
				}

				if(fabs(a[ii1][ii2]) < (10.0*DBL_MIN)){
					a[ii1][ii2] = 0.0;
					continue;
				}

				if( (nearly_eq(fabs(a[ii1][ii1]) + g, fabs(a[ii1][ii1])))
				  &&(nearly_eq(fabs(a[ii2][ii2]) + g, fabs(a[ii2][ii2]))) ){
					a[ii1][ii2] = 0.0;
					continue;
				}

				h = a[ii2][ii2] - a[ii1][ii1];
				if( nearly_eq(fabs(h) + g, fabs(h)) ){
					t = (a[ii1][ii2])/h;
				}else{
					double theta = 0.5*h/(a[ii1][ii2]);

					t = 1.0/(fabs(theta)+sqrt(1.0+theta*theta));
					if(theta < 0.0) t = -t;
				}
				c = 1.0/sqrt(1 + t*t);
				s = t*c;
				tau = s/(1.0 + c);
				t_a = t * a[ii1][ii2];
				add_d[ii1] -= t_a;
				add_d[ii2] += t_a;
				a[ii1][ii1] -= t_a;
				a[ii2][ii2] += t_a;
				a[ii1][ii2] = 0.0;

				for(ii3=0; ii3<ii1; ii3++){
					double tmp_a1 = a[ii3][ii1];
					double tmp_a2 = a[ii3][ii2];

					a[ii3][ii1] = tmp_a1 - s*(tmp_a2 + tmp_a1*tau);
					a[ii3][ii2] = tmp_a2 + s*(tmp_a1 - tmp_a2*tau);
				}

				for(ii3=ii1+1; ii3<ii2; ii3++){
					double tmp_a1 = a[ii1][ii3];
					double tmp_a2 = a[ii3][ii2];

					a[ii1][ii3] = tmp_a1 - s*(tmp_a2 + tmp_a1*tau);
					a[ii3][ii2] = tmp_a2 + s*(tmp_a1 - tmp_a2*tau);
				}

				for(ii3=ii2+1; ii3<n; ii3++){
					double tmp_a1 = a[ii1][ii3];
					double tmp_a2 = a[ii2][ii3];

					a[ii1][ii3] = tmp_a1 - s*(tmp_a2 + tmp_a1*tau);
					a[ii2][ii3] = tmp_a2 + s*(tmp_a1 - tmp_a2*tau);
				}

				for(ii3=0; ii3<n; ii3++){
					double tmp_a1 = e_vects[ii3][ii1];
					double tmp_a2 = e_vects[ii3][ii2];

					e_vects[ii3][ii1] = tmp_a1 - s*(tmp_a2 + tmp_a1*tau);
					e_vects[ii3][ii2] = tmp_a2 + s*(tmp_a1 - tmp_a2*tau);
				}
			}
		}

		/* 対角要素の加算精度を向上 */
		for(ii1=0; ii1<n; ii1++){
			e_vals[ii1] += add_d[ii1];
			a[ii1][ii1] = e_vals[ii1];
			add_d[ii1] = 0.0;
		}
	}

	for(ii1=0; ii1<n; ii1++){
		int max_index = ii1;

		for(ii2=ii1+1; ii2<n; ii2++){
			if(fabs(e_vals[ii2]) > fabs(e_vals[max_index])) max_index = ii2;
		}

		if(max_index != ii1){
			double tmp;

			tmp = e_vals[ii1];
			e_vals[ii1] = e_vals[max_index];
			e_vals[max_index] = tmp;

			for(ii2=0; ii2<n; ii2++){
				tmp = e_vects[ii2][ii1];
				e_vects[ii2][ii1] = e_vects[ii2][max_index];
				e_vects[ii2][max_index] = tmp;
			}
		}
	}

	free((void *)add_d);

	return(NORMAL_RC);
}


RC double2int_pos(double dv, int *iv)
{
	if(dv < -ADM_ABS_ERROR) return(CAL_ERROR_RC);
	if(dv < ADM_ABS_ERROR){
		*iv = 0;
		return(NORMAL_RC);
	}

	*iv = (int)(dv + 0.4999);
	if( !nearly_eq((double)(*iv), dv) ) return(CAL_ERROR_RC);

	return(NORMAL_RC);
}


long double factorial(int i)
{
	int ii1;
	long double ret = 1.0L;

	if(i < 0) return(0.0L);
	for(ii1=1; ii1<=i; ii1++){
		ret *= (long double)ii1;
	}

	return(ret);
}


long double permutation(int n, int r)
{
	int ii1;
	long double ret = 1.0L;

	if( (n < 0) || (r < 0) || (r > n) ) return(0.0L);
	for(ii1=n-r+1; ii1<=n; ii1++){
		ret *= (long double)ii1;
	}

	return(ret);
}


long double combination(int n, int r)
{
	long double ret;

	if( (n < 0) || (r < 0) || (r > n) ) return(0.0L);
	if((n - r) < r) r = n - r;
	ret = permutation(n, r)/factorial(r);

	return(ret);
}


/* A[n][n] を直交行列 Q[n][n] と 上(右)三角R[n][n] に分解 */
RC qr_factor(int n, double **A, double **Q, double **R)
{
	int ii1, ii2, ii3, ii4;
	double tmp_product, max_product;
	double min_norm;
	int *zero_flags;
	double *vect_sum;
	int max_index;

	RC_NULL_CHK( zero_flags = (int *)malloc(sizeof(int)*n) );
	for(ii1=0; ii1<n; ii1++) zero_flags[ii1] = 0;
	RC_NULL_CHK( vect_sum = (double *)malloc(sizeof(double)*n) );

	init_matrix(n, n, R);

	/* 列同士の内積を高速化するため A^t => Q */
	/* 以後，A の列参照は Q の行参照で代用   */
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

	/* zero_flags の立っている R の行をゼロにする */
	for(ii1=0; ii1<n; ii1++){
		if(zero_flags[ii1] == 0) continue;
		for(ii2=ii1; ii2<n; ii2++){
			R[ii1][ii2] = 0.0;
		}
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
	transpose_overwrite_matrix(n, Q);

	free(zero_flags);

	return(NORMAL_RC);
}


RC qr_subst(int n, double **Q, double **R, const double *b, double *x)
{
	int ii1, ii2;

	/* Q^T * b => x */
	for(ii1=0; ii1<n; ii1++){
		x[ii1] = 0.0;
	}
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			x[ii2] += Q[ii1][ii2] * b[ii1];
		}
	}

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


RC qr_solve(double **A, int n, const double *b, double *x)
{
	double **Q;
	double **R;

	RC_TRY( allocate2D(n, n, &Q) );
	RC_TRY( allocate2D(n, n, &R) );

	RC_TRY( qr_factor(n, A, Q, R) );
	RC_TRY( qr_subst(n, Q, R, b, x) );

	RC_TRY( free2D(n, n, Q) );
	RC_TRY( free2D(n, n, R) );

	return(NORMAL_RC);
}


void my_sort_int(int n, int array[])
{
	int ii1, ii2;
	int tmp;

	for(ii1=0; ii1<n; ii1++){
		int min_index = ii1;
		for(ii2=ii1+1; ii2<n; ii2++){
			if(array[ii2] < array[min_index]) min_index = ii2;
		}
		tmp = array[min_index];
		array[min_index] = array[ii1];
		array[ii1] = tmp;
	}
}

