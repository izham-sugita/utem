/*********************************************************************
 * matrix.c
 *
 * Copyright (C) 2000 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Syoji ITO> <Masanobu SUYAMA> <Ryu SASAOKA>
 *  <Takaaki NAGATANI> <Yasuo ASAGA>
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: matrix.c 1066 2016-05-19 09:53:33Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"

static void mul_matrix_large(int l, int m, int n,
                             double **A, double **B, double **C);
static void mul_matrix_AtBA3(int m, int n, double **A, double **B, double **C);
static double determinant2(double **matrix);
static double determinant3(double **matrix);


RC
allocate1D (int size, double **array)
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (double *)mm_alloc(size * sizeof(double)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1] = 0.0;
	}

	return(NORMAL_RC);
}


RC
allocate2D (int record, int column, double ***array)
{
	int ii1;
	double *ptr;

	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (double **)mm_alloc(record * sizeof(double *)) );
	RC_TRY( allocate1D(record*column, &ptr) );

	for(ii1=0; ii1<record; ii1++){
		(*array)[ii1] = &(ptr[ii1*column]);
	}

	return(NORMAL_RC);
}


RC
allocate3D (int record1, int record2, int record3, double ****array)
{
	int ii1;
	double **ptr;

	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (double ***)mm_alloc(record1 * sizeof(double **)) );
	RC_TRY( allocate2D(record1*record2, record3, &ptr) );

	for(ii1=0;ii1<record1;ii1++){
		(*array)[ii1] = &(ptr[ii1*record2]);
	}

	return(NORMAL_RC);
}


RC
allocate1I (int size, int **array)
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (int *)mm_alloc(size * sizeof(int)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1] = 0;
	}

	return(NORMAL_RC);
}


RC
allocate2I (int record, int column, int ***array)
{
	int ii1;
	int *ptr;

	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (int **)mm_alloc(record * sizeof(int *)) );
	RC_TRY( allocate1I(record*column, &ptr) );

	for(ii1=0; ii1<record; ii1++){
		(*array)[ii1] = &(ptr[ii1*column]);
	}

	return(NORMAL_RC);
}


RC
allocate3I (int record1, int record2, int record3, int ****array)
{
	int ii1;
	int **ptr;

	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (int ***)mm_alloc(record1 * sizeof(int **)) );
	RC_TRY( allocate2I(record1*record2, record3, &ptr) );

	for(ii1=0;ii1<record1;ii1++){
		(*array)[ii1] = &(ptr[ii1*record2]);
	}

	return(NORMAL_RC);
}



RC
free1D (int size, double **array)
{
	if(size <= 0) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free2D (int record, int column, double ***array)
{
	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free3D (int record1, int record2, int record3, double ****array)
{
	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((**array)[0]) );
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free1I (int size, int **array)
{
	if(size <= 0) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free2I (int record, int column, int ***array)
{
	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC
free3I (int record1, int record2, int record3, int ****array)
{
	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((**array)[0]) );
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


/* array[record*column] を matrix[record][column] に変換 */
RC
array2matrix (int record, int column, double array[], double *matrix[])
{
	int ii1;

	for(ii1=0; ii1<record; ii1++){
		matrix[ii1] = &(array[ii1*column]);
	}

	return(NORMAL_RC);
}


/* wa * A[m][n] => B[m][n] */
/* A == B でも構わない */
void
mul_scalar_matrix (int m, int n, double wa, double **A, double **B)
{
	int ii1, ii2;

	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n; ii2++){
			B[ii1][ii2] = wa * A[ii1][ii2];
		}
	}
}


static void
mul_matrix_large (int l, int m, int n, double **A, double **B, double **C)
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
void
mul_matrix (int l, int m, int n, double **A, double **B, double **C)
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
void
mul_matrix_vect (int m, int n, double **A, const double *X, double *Y)
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
void
mul_trans_matrix_vect (int m, int n, double **A, const double *X, double *Y)
{
	int ii1, ii2;

	for(ii1=0; ii1<n; ii1++) Y[ii1] = 0.0;
	for(ii1=0; ii1<m; ii1++){
		for(ii2=0; ii2<n; ii2++){
			Y[ii2] += A[ii1][ii2] * X[ii1];
		}
	}
}


/* X^T[m] * A[m][n] => Y[n] */
void
mul_vect_matrix (int m, int n, const double *X, double **A, double *Y)
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
double
determinant (int dim, double **matrix)
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
static double
determinant2 (double **matrix)
{
	return(matrix[0][0] * matrix[1][1]
	     - matrix[0][1] * matrix[1][0]);
}


/* 3 x 3 行列の行列式 */
static double
determinant3 (double **matrix)
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
RC
inverse_matrix (int dim, double **A, double **inverseA)
{
	double div;

	if(dim < 1) return(ARG_ERROR_RC);

	if(dim == 1){
		div = keep_away_zero(A[0][0]);
		inverseA[0][0] = 1.0/div;
	}else if(dim == 2){
		div = keep_away_zero(determinant2(A));
		inverseA[0][0] = A[1][1] / div;
		inverseA[0][1] = - A[0][1] / div;
		inverseA[1][0] = - A[1][0] / div;
		inverseA[1][1] = A[0][0] / div;
	}else if(dim == 3){
		div = keep_away_zero(determinant3(A));
		inverseA[0][0] = (- A[1][2] * A[2][1] + A[1][1] * A[2][2])/div;
		inverseA[0][1] = (  A[0][2] * A[2][1] - A[0][1] * A[2][2])/div;
		inverseA[0][2] = (- A[0][2] * A[1][1] + A[0][1] * A[1][2])/div;
		inverseA[1][0] = (  A[1][2] * A[2][0] - A[1][0] * A[2][2])/div;
		inverseA[1][1] = (- A[0][2] * A[2][0] + A[0][0] * A[2][2])/div;
		inverseA[1][2] = (  A[0][2] * A[1][0] - A[0][0] * A[1][2])/div;
		inverseA[2][0] = (- A[1][1] * A[2][0] + A[1][0] * A[2][1])/div;
		inverseA[2][1] = (  A[0][1] * A[2][0] - A[0][0] * A[2][1])/div;
		inverseA[2][2] = (- A[0][1] * A[1][0] + A[0][0] * A[1][1])/div;
	}else{
		int ii1, ii2;
		double dum;
		double *col;
		double **LU;
		int *index;

		col = (double *)mm_alloc(dim * sizeof(double));
		index = (int *)mm_alloc(dim * sizeof(int));
		if((col == NULL)||(index == NULL)){
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

		RC_TRY( mm_free(col) );
		RC_TRY( mm_free(index) );
		RC_TRY( free2D(dim, dim, &LU) );
	}

	return(NORMAL_RC);
}


/* wa * A[m][n] + wb * B[m][n] => C[m][n] */
/* A == C or B == C でも構わない          */
RC
add_matrix (int m, int n, double wa, double **A, double wb, double **B,
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

RC copy_vect (int size, double *vect1, double *vect2)
{
	RC_NEG_ZERO_CHK( size );
	RC_NULL_CHK( vect1 );
	RC_NULL_CHK( vect2 );

	int ii1;
	for(ii1=0; ii1<size; ii1++){
		vect2[ii1] = vect1[ii1];
	}

	return(NORMAL_RC);
}

RC copy_matrix (int m, int n, double **A, double **B)
{
	RC_NEG_ZERO_CHK( m );
	RC_NEG_ZERO_CHK( n );
	RC_NULL_CHK( A );
	RC_NULL_CHK( B );
	
	int ii1;
	for(ii1=0; ii1<m; ii1++){
		RC_TRY( copy_vect(n, A[ii1], B[ii1]) );
	}

	return(NORMAL_RC);
}


/* 行列 A (row * col) を転置して B に代入 */
void
transpose_matrix (int row, int col, double **A, double **B)
{
	int ii1, ii2;

	for(ii1=0; ii1<col; ii1++){
		for(ii2=0; ii2<row; ii2++){
			B[ii1][ii2] = A[ii2][ii1];
		}
	}
}

/* 正方行列 A[n][n] を転置して上書き */
void
transpose_overwrite_matrix (int n, double **A)
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


/* ベクトル v[n] を 0.0 で初期化する．
 */
void
init_vect (int n, double v[])
{
	int ii1;

	for(ii1=0; ii1<n; ii1++){
		v[ii1] = 0.0;
	}
}


/* 行列 A[row][col] を 0.0 で初期化する． 
 */
void
init_matrix (int row, int col, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<row; ii1++){
		for(ii2=0; ii2<col; ii2++){
			A[ii1][ii2] = 0.0;
		}
	}
}

/* ベクトル v[n] を単位ベクトルにする
 */
RC unit_vect (int n, double *v)
{
	int ii1;
	double norm = norm_vect(n, v);
	if(norm < ABS_ERROR) return(CAL_ERROR_RC);

	for(ii1=0; ii1<n; ii1++){
		v[ii1] /= norm;
	}

	return(NORMAL_RC);
}


/* 行列 A (n * n) を単位行列として初期化 */
void
unit_matrix (int n, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			A[ii1][ii2] = 0.0;
		}
		A[ii1][ii1] = 1.0;
	}
}


/* v1[m] と v2[n] のテンソル積 => A[m][n] */
RC
tensor_product (int m, int n, const double *v1, const double *v2, double **A)
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
double
inner_product (int n, const double *v1, const double *v2)
{
	long double ret = 0.0;
	int ii1;

	for(ii1=0; ii1<n; ii1++){
		ret += v1[ii1]*v2[ii1];
	}

	return((double)ret);
}


/* v1[n] と v2[n] の外積を計算する．
 * 3次元のみ．
 */
RC outer_product (const double *v1, const double *v2, double *ret)
{
	RC_NULL_CHK( ret );

	ret[0] = v1[1]*v2[2] - v1[2]*v2[1];
	ret[1] = v1[2]*v2[0] - v1[0]*v2[2];
	ret[2] = v1[0]*v2[1] - v1[1]*v2[0];

	return(NORMAL_RC);
}


/* v1[n] と v2[n] の内積(桁落ち対策, long double 版) */
long double
inner_product_ac(int n, const long double v1[], const long double v2[])
{
	long double ret, r;
	long double vtmp;
	int ii1;

	r = 0.0;
	ret = 0.0;
	for(ii1=0; ii1<n-5; ii1+=6){
		ERROR_SUM( ret, v1[ii1  ]*v2[ii1  ] + v1[ii1+1]*v2[ii1+1]
		              + v1[ii1+2]*v2[ii1+2] + v1[ii1+3]*v2[ii1+3]
		              + v1[ii1+4]*v2[ii1+4] + v1[ii1+5]*v2[ii1+5], r, vtmp );
	}
	for(; ii1<n; ii1++){
		ERROR_SUM( ret, v1[ii1]*v2[ii1], r, vtmp );
	}

	return(ret);
}

/* ベクトルのL2ノルム */
double norm_vect(int size, double *vect)
{
	return(sqrt(inner_product(size, vect, vect)));
}

/* A[l][m]^T * B[l][n] => C[m][n] */
void
mul_matrix_AtB (int l, int m, int n, double **A, double **B, double **C)
{
	int ii1, ii2, ii3;

	if((l <= 0)||(m <= 0)||(n <= 0)) return;
	for(ii1=0; ii1<m; ii1++){
		double *C_ii1 = C[ii1];
		double A_0_ii1 = A[0][ii1];
		double *B_0 = B[0];

		for(ii3=0; ii3<n; ii3++){
			C_ii1[ii3] = A_0_ii1 * B_0[ii3];
		}

		for(ii2=1; ii2<l; ii2++){
			double A_ii2_ii1 = A[ii2][ii1];
			double *B_ii2 = B[ii2];

			for(ii3=0; ii3<n; ii3++){
				C_ii1[ii3] += A_ii2_ii1 * B_ii2[ii3];
			}
		}
	}
}

/* A[m][n]^T * B[m][m] * A[m][n] = C[n][n] */
void
mul_matrix_AtBA (int m, int n, double **A, double **B, double **C)
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
static void
mul_matrix_AtBA3 (int m, int n, double **A, double **B, double **C)
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


/* print v[n] */
void
print_vect (FILE *fp, int n, const double v[])
{
	int ii1;

	for(ii1=0; ii1<n; ii1++){
		fprintf(fp, " %11.3e", v[ii1]);
	}
	fprintf(fp, "\n");
}


/* print A[m][n] */
void
print_matrix (FILE *fp, int m, int n, double **A)
{
	int ii1;

	for(ii1=0; ii1<m; ii1++){
		print_vect(fp, n, A[ii1]);
	}
}


/*
 * a[n] = A[n][n]*b[n]
 * a[] の内，flag_a[n] = 1 の変数を既知，= 0 の変数を未知
 * b[] の内，flag_b[n] = 1 の変数を既知，= 0 の変数を未知
 * a[], b[] の既知変数の数は合計 n とする．
 * この時，y[n] = A[n][n]*x[n] となるように A[][] を書き換える
 * y[n] : a[] の未知変数，b[] の未知変数を順に並べたベクトル
 * x[n] : a[] の既知変数，b[] の既知変数を順に並べたベクトル 
 * tmp_matrix[7][n][n]
 */
RC
flag_convert_matrix (int n, double **A, int flag_a[], int flag_b[],
                     double ***tmp_matrix)
{
	int ii1, ii2;
	int a_count, b_count;
	double **A00 = tmp_matrix[0];
	double **A01 = tmp_matrix[1];
	double **A10 = tmp_matrix[2];
	double **A11 = tmp_matrix[3];
	double **A10inv = tmp_matrix[4];
	double **Atmp0 = tmp_matrix[5];
	double **Atmp1 = tmp_matrix[6];
	int wp1, wp2;

	if(n <= 0) return(ARG_ERROR_RC);
	a_count = b_count = 0;
	for(ii1=0; ii1<n; ii1++){
		if(flag_a[ii1] == 0) a_count++;
		if(flag_b[ii1] == 0) b_count++;
	}
	if((a_count + b_count) != n) return(ARG_ERROR_RC);
	if(b_count == 0) return(NORMAL_RC);

	/* A00 */
	wp1 = wp2 = 0;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			if((flag_a[ii1] == 0)&&(flag_b[ii2] == 0)){
				A00[wp1][wp2] = A[ii1][ii2];
				wp2++;
			}
		}
		if(wp2 > 0) wp1++;
		wp2 = 0;
	}

	/* A01 */
	wp1 = wp2 = 0;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			if((flag_a[ii1] == 0)&&(flag_b[ii2] == 1)){
				A01[wp1][wp2] = A[ii1][ii2];
				wp2++;
			}
		}
		if(wp2 > 0) wp1++;
		wp2 = 0;
	}

	/* A10 */
	wp1 = wp2 = 0;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			if((flag_a[ii1] == 1)&&(flag_b[ii2] == 0)){
				A10[wp1][wp2] = A[ii1][ii2];
				wp2++;
			}
		}
		if(wp2 > 0) wp1++;
		wp2 = 0;
	}

	/* A11 */
	wp1 = wp2 = 0;
	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<n; ii2++){
			if((flag_a[ii1] == 1)&&(flag_b[ii2] == 1)){
				A11[wp1][wp2] = A[ii1][ii2];
				wp2++;
			}
		}
		if(wp2 > 0) wp1++;
		wp2 = 0;
	}

	RC_TRY( inverse_matrix(b_count, A10, A10inv) );
	if(a_count == 0){
		RC_TRY( copy_matrix(n, n, A10inv, A) );
		return(NORMAL_RC);
	}

	mul_matrix(a_count, b_count, b_count, A00, A10inv, Atmp0);
	for(ii1=0; ii1<a_count; ii1++){
		for(ii2=0; ii2<b_count; ii2++){
			A[ii1][ii2] = Atmp0[ii1][ii2];
		}
	}

	mul_matrix(b_count, b_count, a_count, A10inv, A11, Atmp0);
	mul_matrix(a_count, b_count, a_count, A00, Atmp0, Atmp1);
	add_matrix(a_count, a_count, 1.0, A01, -1.0, Atmp1, Atmp0);
	for(ii1=0; ii1<a_count; ii1++){
		for(ii2=0; ii2<a_count; ii2++){
			A[ii1][b_count + ii2] = Atmp0[ii1][ii2];
		}
	}

	for(ii1=0; ii1<b_count; ii1++){
		for(ii2=0; ii2<b_count; ii2++){
			A[a_count + ii1][ii2] = A10inv[ii1][ii2];
		}
	}

	mul_matrix(b_count, b_count, a_count, A10inv, A11, Atmp0);
	for(ii1=0; ii1<b_count; ii1++){
		for(ii2=0; ii2<a_count; ii2++){
			A[a_count + ii1][b_count + ii2] = -Atmp0[ii1][ii2];
		}
	}

	return(NORMAL_RC);
}


/* A[n][n] の対称成分で A[n][n] を上書き */
RC
symmetrize_matrix (int n, double **A)
{
	int ii1, ii2;

	for(ii1=0; ii1<(n-1); ii1++){
		for(ii2=0; ii2<ii1; ii2++){
			A[ii1][ii2] = A[ii2][ii1] = (A[ii1][ii2] + A[ii2][ii1])/2.0;
		}
	}

	return(NORMAL_RC);
}


/* A[n][n] の j 行，j 列の w? 倍を i? 行，i? 列に加算し，*/
/* j 行，j 列をゼロにする                                */
/* y != NULL なら右辺ベクトルとして同時に処理            */
RC
row_col_add_zero3_matrix (int n, double **A, double y[], int i1, double w1,
                          int i2, double w2, int i3, double w3, int j)
{
	int ii1;

	if( (n <= 0)||(j < 0)||(j >= n) ) return(ARG_ERROR_RC);
	if(i1 >= n) return(ARG_ERROR_RC);
	if(i2 >= n) return(ARG_ERROR_RC);
	if(i3 >= n) return(ARG_ERROR_RC);

	if(i1 >= 0){
		for(ii1=0; ii1<n; ii1++) A[i1][ii1] += w1*A[j][ii1];
		if(y != NULL) y[i1] += w1*y[j];
	}
	if(i2 >= 0){
		for(ii1=0; ii1<n; ii1++) A[i2][ii1] += w2*A[j][ii1];
		if(y != NULL) y[i2] += w2*y[j];
	}
	if(i3 >= 0){
		for(ii1=0; ii1<n; ii1++) A[i3][ii1] += w3*A[j][ii1];
		if(y != NULL) y[i3] += w3*y[j];
	}
	for(ii1=0; ii1<n; ii1++) A[j][ii1] = 0.0;
	if(y != NULL) y[j] = 0.0;

	if(i1 >= 0){
		for(ii1=0; ii1<n; ii1++) A[ii1][i1] += w1*A[ii1][j];
	}
	if(i2 >= 0){
		for(ii1=0; ii1<n; ii1++) A[ii1][i2] += w2*A[ii1][j];
	}
	if(i3 >= 0){
		for(ii1=0; ii1<n; ii1++) A[ii1][i3] += w3*A[ii1][j];
	}
	for(ii1=0; ii1<n; ii1++) A[ii1][j] = 0.0;

	return(NORMAL_RC);
}


RC
row_col_add_zero_matrix (int n, double **A, double y[], int i_size,
                         const int i[], const double w[], int j)
{
	int ii1, ii2;

	if( (n <= 0)||(j < 0)||(j >= n) ) return(ARG_ERROR_RC);
	for(ii1=0; ii1<i_size; ii1++){
		if(i[ii1] >= n) return(ARG_ERROR_RC);
	}

	for(ii1=0; ii1<i_size; ii1++){
		if(i[ii1] >= 0){
			for(ii2=0; ii2<n; ii2++) A[i[ii1]][ii2] += w[ii1]*A[j][ii2];
			if(y != NULL) y[i[ii1]] += w[ii1]*y[j];
		}
	}
	for(ii1=0; ii1<n; ii1++) A[j][ii1] = 0.0;
	if(y != NULL) y[j] = 0.0;

	for(ii1=0; ii1<i_size; ii1++){
		if(i[ii1] >= 0){
			for(ii2=0; ii2<n; ii2++) A[ii2][i[ii1]] += w[ii1]*A[ii2][j];
		}
	}
	for(ii1=0; ii1<n; ii1++) A[ii1][j] = 0.0;

	return(NORMAL_RC);
}


/* A[n][n] のマトリックスを静的縮約により i 行, i 列を消去 */
/* y != NULL なら右辺ベクトルとして同時に処理              */
RC
static_condens_matrix (int n, double **A, double y[], int i)
{
	int ii1, ii2;
	long double inv;

	if( (n <= 0)||(i < 0)||(i >= n) ) return(ARG_ERROR_RC);

	inv = 1.0/(A[i][i] + SIGN(A[i][i])*ABS_TOL);
	for(ii1=0; ii1<n; ii1++){
		if(ii1 == i) continue;
		for(ii2=0; ii2<i; ii2++){
			A[ii1][ii2] -= inv*A[ii1][i]*A[i][ii2];
		}
		for(ii2=i+1; ii2<n; ii2++){
			A[ii1][ii2] -= inv*A[ii1][i]*A[i][ii2];
		}
		if(y != NULL) y[ii1] -= inv*A[ii1][i]*y[i];
	}

	for(ii1=0; ii1<n; ii1++){
		A[ii1][i] = A[i][ii1] = 0.0;
	}
	if(y != NULL) y[i] = 0.0;

	return(NORMAL_RC);
}

#define J_TRANS_MAX 100
/* Jacobi 変換による対称行列a[n][n]の固有値e_vals[n]、*/
/* 固有ベクトルe_vects[n][n]の計算                    */
/* a[][] は対角行列に書き換えられえる                 */
RC
eigen_jacobi_trans (int n, double **a, double *e_vals, double **e_vects)
{
	int ii1, ii2, ii3, trans_count, trans_max;
	double sum, sum_d, thresh;
	double *add_d; /* 対角要素に加算する数値 */

	if(n < 1) return(ARG_ERROR_RC);
	trans_max = J_TRANS_MAX + n;

	unit_matrix(n, e_vects);

	add_d = (double *)mm_alloc(n * sizeof(double));
	if(add_d == NULL) return(ALLOC_ERROR_RC);
	for(ii1=0; ii1<n; ii1++){
		e_vals[ii1] = a[ii1][ii1];
		add_d[ii1] = 0.0;
	}

	trans_count = 0;
	while(1){
		if(trans_count > trans_max) return(CAL_ERROR_RC);
		trans_count++;

		sum = 0.0;
		sum_d = 0.0;
		for(ii1=0; ii1<n; ii1++){
			for(ii2=ii1+1; ii2<n; ii2++){
				sum += fabs(a[ii1][ii2]);
			}
			sum_d += fabs(a[ii1][ii1]);
		}

		if(sum <= (10.0*DBL_MIN)) break; /* 絶対収束 */
		
		if((sum/(double)n) <= ((10.0*DBL_EPSILON)*sum_d)){ /* 収束判定 */
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

	RC_TRY( mm_free((void *)add_d) );

	return(NORMAL_RC);
}


/* b[n] = A[n][n]*x[n] */
/* x[0] ～ x[nknown-1] : 既知 */
/* x[nknown] ～ x[n-1] : 未知 */
/* 本関数呼び出し時に b[] の既知部分は x[] の未知部分に代入しておく */
/* x[] の未知部分を計算して上書きする */
RC
simple_subst_solve_matrix (int n, double **A, double x[], int nknown,
                           double **tmp_matrix)
{
	int ii1, ii2;


	if(n <= 0) return(ARG_ERROR_RC);
	if( (nknown < 0)||(nknown > n) ) return(ARG_ERROR_RC);
	if(nknown == n) return(NORMAL_RC);

	/* x[未] -= A[未][既] * x[既] */
	for(ii1=nknown; ii1<n; ii1++){
		for(ii2=0; ii2<nknown; ii2++){
			x[ii1] -= A[ii1][ii2] * x[ii2];
		}
	}

	/* tmp_matrix = A[未][未] */
	for(ii1=nknown; ii1<n; ii1++){
		for(ii2=nknown; ii2<n; ii2++){
			tmp_matrix[ii1 - nknown][ii2 - nknown] = A[ii1][ii2];
		}
	}

	/* x[未] = (A[未][未])^-1 * x[未] */
	RC_TRY( lu_solve(tmp_matrix, n - nknown, &(x[nknown]), &(x[nknown])) );

	return(NORMAL_RC);
}


