/*********************************************************************
 * complex_utl.c
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA> <Kohei SHINTANI> <Yuri NAKAMURA> 
 *  <Takuya HAYASHI>
 *
 *  Refer to log documents for details.
 *********************************************************************/

/* $Id: complex_utl.c 1090 2016-12-01 03:23:36Z hayashi $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "memory_manager.h"
#include "math_utl.h"
#include "complex_utl.h"


/* 極座標系表示 (r, theta) からDescartes座標系表示 (x, y) を計算する．
 */
COMPLEX polar2right_complex (double r, double theta)
{
	return(r*(cos(theta) + sin(theta)*COMPLEX_UNIT));
}


/* 複素数型の1次元配列を動的確保する．
 */
RC allocate1C(int size, COMPLEX **array)
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (COMPLEX *)mm_alloc(size * sizeof(COMPLEX)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1] = 0.0;
	}

	return(NORMAL_RC);
}

/* 複素数型の2次元配列を動的確保する．
 */
RC allocate2C(int record, int column, COMPLEX ***array)
{
	int ii1;
	COMPLEX *ptr;

	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array)
	              = (COMPLEX **)mm_alloc(record * sizeof(COMPLEX *)) );
	RC_TRY( allocate1C(record*column, &ptr) );

	for(ii1=0; ii1<record; ii1++){
		(*array)[ii1] = &(ptr[ii1*column]);
	}

	return(NORMAL_RC);
}

/* 複素数型の3次元配列を動的確保する．
 */
RC allocate3C(int record1, int record2, int record3, COMPLEX ****array)
{
	int ii1;
	COMPLEX **ptr;

	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array)
	              = (COMPLEX ***)mm_alloc(record1 * sizeof(COMPLEX **)) );
	RC_TRY( allocate2C(record1*record2, record3, &ptr) );

	for(ii1=0;ii1<record1;ii1++){
		(*array)[ii1] = &(ptr[ii1*record2]);
	}

	return(NORMAL_RC);
}

/* 動的確保した複素数型の1次元配列を解放する．
 */
RC free1C(int size, COMPLEX **array)
{
	if(size <= 0) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}

/* 動的確保した複素数型の1次元配列を解放する．
 */
RC free2C(int record, int column, COMPLEX ***array)
{
	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}

/* 動的確保した複素数型の1次元配列を解放する．
 */
RC free3C (int record1, int record2, int record3, COMPLEX ****array)
{
	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((**array)[0]) );
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


/* 要素数 size の1次元配列 vect1 を vect2 にコピーする．
 * vect2 は動的確保してから渡す．
 */
RC copy_vect_complex (int size, const COMPLEX *vect1, COMPLEX *vect2)
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

/* 要素数 size の1次元配列 vect のすべての要素に 0.0 を代入する．
 */
RC init_vect_complex (int size, COMPLEX *vect)
{
	RC_NEG_ZERO_CHK( size );
	RC_NULL_CHK( vect );

	int ii1;
	for(ii1=0; ii1<size; ii1++){
		vect[ii1] = 0.0;
	}

	return(NORMAL_RC);
}

/* 行数 row，列数 col の2次元配列 matrix のすべての要素に 0.0 を代入する．
 */
RC init_matrix_complex (int row, int col, COMPLEX **matrix)
{
	RC_NEG_ZERO_CHK( row );
	RC_NEG_ZERO_CHK( col );
	RC_NULL_CHK( matrix );

	int ii1;
	for(ii1=0; ii1<row; ii1++){
		RC_TRY( init_vect_complex(col, matrix[ii1]) );
	}

	return(NORMAL_RC);
}

/* 要素数 size の1次元配列 vect を fp に出力する．
 */
RC print_vect_complex (FILE *fp, int size, const COMPLEX *vect)
{
	RC_NULL_CHK( fp );
	RC_NEG_ZERO_CHK( size );
	RC_NULL_CHK( vect );

	int ii1;
	for(ii1=0; ii1<size; ii1++){
		if(ii1 > 0) fprintf(fp, ", ");
		fprintf(fp, "%5.3e+%5.3ei", creal(vect[ii1]), cimag(vect[ii1]));
	}
	fprintf(fp, "\n");

	return(NORMAL_RC);
}

/* 行数 row，列数 col の2次元配列 matrix を fp に出力する．
 */
RC print_matrix_complex (FILE *fp, int row, int col, const COMPLEX **matrix)
{
	RC_NULL_CHK( fp );
	RC_NEG_ZERO_CHK( row );
	RC_NEG_ZERO_CHK( col );
	RC_NULL_CHK( matrix );

	int ii1;
	for(ii1=0; ii1<row; ii1++){
		print_vect_complex(fp, col, matrix[ii1]);
	}

	return(NORMAL_RC);
}

/* 複素行列 matrix と複素ベクトル vect の積を計算する．
 *   ret = matrix[row][col] * vect[col] 
 */
void mul_matrix_vect_complex (int row, int col, COMPLEX **matrix, 
		COMPLEX *vect, COMPLEX *ret)
{
	int ii1, ii2;
	long double _Complex sum;
	for(ii1=0; ii1<row; ii1++){
		sum = 0.0;

		for(ii2=0; ii2<col; ii2++){
			sum += matrix[ii1][ii2]*vect[ii2];
		}
		ret[ii1] = (COMPLEX)sum;
	}
}

/* 複素ベクトル vect の Euclidノルムを計算する． 
 */
double norm_vect_complex (int size, COMPLEX *vect)
{
	return((double)csqrt(hermite_product(size, vect, vect)));
}

/* 2つの複素ベクトル vect1 と vect2 のHermite積を計算する．
 */
COMPLEX hermite_product (int size, COMPLEX *vect1, COMPLEX *vect2)
{
	int ii1;
	COMPLEX ret = 0.0;
	for(ii1=0; ii1<size; ii1++){
		ret += conj(vect1[ii1])*vect2[ii1];
	}

	return(ret);
}


