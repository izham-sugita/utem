/*********************************************************************
 * complex_utl.c
 *
 * Copyright (C) 2009 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA> <Kohei SHINTANI> <Yuri NAKAMURA>
 *
 *  Refer to log documents for details.
 *********************************************************************/

/* $Id: complex_utl.c 931 2009-12-14 05:32:09Z aoyama $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rc.h"
#include "math_utl.h"
#include "memory_manager.h"


#define SUB_MUL_COMPLEX(a, b, c) \
a.re = a.re - (b.re * c.re - b.im * c.im); \
a.im = a.im - (b.re * c.im + b.im * c.re);


COMPLEX polar2right_complex(double r, double theta)
{
	COMPLEX ret;

	ret.re = r*cos(theta);
	ret.im = r*sin(theta);

	return(ret);
}


/* a + b */
COMPLEX add_complex(COMPLEX a, COMPLEX b)
{
	COMPLEX ret;

	ret.re = a.re + b.re;
	ret.im = a.im + b.im;

	return(ret);
}


/* a - b */
COMPLEX sub_complex(COMPLEX a, COMPLEX b)
{
	COMPLEX ret;

	ret.re = a.re - b.re;
	ret.im = a.im - b.im;

	return(ret);
}


/* a * b */
COMPLEX mul_complex(COMPLEX a, COMPLEX b)
{
	COMPLEX ret;

	ret.re = a.re * b.re - a.im * b.im;
	ret.im = a.re * b.im + a.im * b.re;

	return(ret);
}


/* a / b */
/* by Iri and Fujino */
COMPLEX div_complex(COMPLEX a, COMPLEX b)
{
	COMPLEX ret;

	b.re = keep_away_zero(b.re);
	b.im = keep_away_zero(b.im);

	if(fabs(b.re) >= fabs(b.im)){
		double div = b.im/b.re;
		double s = b.re + b.im*div;

		ret.re = (a.re + a.im*div)/s;
		ret.im = (-a.re*div + a.im)/s;
	}else{
		double div = b.re/b.im;
		double s = b.re*div + b.im;

		ret.re = (a.re*div + a.im)/s;
		ret.im = (-a.re + a.im*div)/s;
	}

	return(ret);
}


/* a* */
COMPLEX conjugate_complex(COMPLEX a)
{
	COMPLEX ret;

	ret.re = a.re;
	ret.im = -a.im;

	return(ret);
}


/* |a| */
/* by Iri and Fujino */
double abs_complex(COMPLEX a)
{
	double ret;
	double abs_re = fabs(a.re);
	double abs_im = fabs(a.im);

	if(abs_im < ABS_TOL){
		ret = abs_re;
	}else if(abs_re < ABS_TOL){
		ret = abs_im;
	}else if(abs_re > abs_im){
		double div = a.im/a.re;
		ret = abs_re * sqrt(1.0 + div*div);
	}else{
		double div = a.re/a.im;
		ret = abs_im * sqrt(1.0 + div*div);
	}

	return(ret);
}


/* |a|^2 */
double sq_abs_complex(COMPLEX a)
{
	return(a.re*a.re + a.im*a.im);
}


RC allocate1C(int size, COMPLEX **array)
{
	int ii1;

	if(size <= 0) return(ARG_ERROR_RC);

	RC_NULL_CHK( (*array) = (COMPLEX *)mm_alloc(size * sizeof(COMPLEX)) );
	for(ii1=0; ii1<size; ii1++){
		(*array)[ii1].re = 0.0;
		(*array)[ii1].im = 0.0;
	}

	return(NORMAL_RC);
}


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


RC free1C(int size, COMPLEX **array)
{
	if(size <= 0) return(ARG_ERROR_RC);

	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC free2C(int record, int column, COMPLEX ***array)
{
	if((record <= 0)||(column <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


RC free3C(int record1, int record2, int record3, COMPLEX ****array)
{
	if((record1 <= 0)||(record2 <= 0)||(record3 <= 0)) return(ARG_ERROR_RC);
	RC_TRY( mm_free((**array)[0]) );
	RC_TRY( mm_free((*array)[0]) );
	RC_TRY( mm_free(*array) );
	*array = NULL;

	return(NORMAL_RC);
}


/* 行列 A[n][n] を LU 分解して A に上書き */
/* index[n] : ピボット選択による行交換の情報 */
RC lu_decomp_complex(COMPLEX **A, int n, int *index)
{
	int ii1, ii2, ii3;
	double *scale_factor;
	double max;
	double tmp;
	COMPLEX ctmp;
	int max_index;

	if( (n <= 0)||(A == NULL)||(index == NULL) ) return(ARG_ERROR_RC);

	/* scale_factor を設定 */
	RC_NULL_CHK( scale_factor = (double *)mm_alloc(n * sizeof(double)) );
	for(ii1=0; ii1<n; ii1++){
		max = 0.0;
		for(ii2=0; ii2<n; ii2++){
			tmp = abs_complex(A[ii1][ii2]);
			if(tmp > max) max = tmp;
		}
		scale_factor[ii1] = 1.0/keep_away_zero(max);
	}

	for(ii1=0; ii1<n; ii1++){
		for(ii2=0; ii2<ii1; ii2++){
			COMPLEX sum = A[ii2][ii1];

			for(ii3=0; ii3<ii2; ii3++){
				SUB_MUL_COMPLEX(sum, A[ii2][ii3], A[ii3][ii1]);
				/* sum -= A[ii2][ii3] * A[ii3][ii1]; */
			}
			A[ii2][ii1] = sum;
		}

		max = 0.0;
		max_index = ii1;
		for(ii2=ii1; ii2<n; ii2++){
			COMPLEX sum = A[ii2][ii1];

			for(ii3=0; ii3<ii1; ii3++){
				SUB_MUL_COMPLEX(sum, A[ii2][ii3], A[ii3][ii1]);
				/* sum -= A[ii2][ii3] * A[ii3][ii1]; */
			}
			A[ii2][ii1] = sum;

			tmp = scale_factor[ii2] * abs_complex(sum);
			if(tmp > max){
				max = tmp;
				max_index = ii2;
			}
		}

		if(max_index != ii1){
			for(ii2=0; ii2<n; ii2++){
				ctmp = A[max_index][ii2];
				A[max_index][ii2] = A[ii1][ii2];
				A[ii1][ii2] = ctmp;
			}
			scale_factor[max_index] = scale_factor[ii1];
		}
		index[ii1] = max_index;

		for(ii2=ii1+1; ii2<n; ii2++){
			A[ii2][ii1] = div_complex(A[ii2][ii1], A[ii1][ii1]);
			/* A[ii2][ii1] /= A[ii1][ii1]; */
		}
	}

	RC_TRY( mm_free(scale_factor) );

	return(NORMAL_RC);
}


/* 前進代入、後退代入を行なう */
/* A[n][n], index[n] : lu_decomp_complex で得た値 */
/* B[n] : 右辺ベクトル -> 解ベクトル */
RC lu_subst_complex(COMPLEX **A, int n, int *index, COMPLEX *B)
{
	int ii1, ii2;
	COMPLEX tmp;

	if((n <= 0)||(A == NULL)||(index == NULL)||(B == NULL)){
		return(ARG_ERROR_RC);
	}

	/* B[] の並び替え、nzero_index の設定 */
	for(ii1=0; ii1<n; ii1++){
		if((index[ii1] < 0)||(index[ii1] >= n)) return(ARG_ERROR_RC);

		tmp = B[index[ii1]];
		B[index[ii1]] = B[ii1];
		B[ii1] = tmp;
	}

	/* 前進代入 */
	for(ii1=1; ii1<n; ii1++){
		for(ii2=0; ii2<ii1; ii2++){
			SUB_MUL_COMPLEX(B[ii1], A[ii1][ii2], B[ii2]);
			/* B[ii1] -= A[ii1][ii2]*B[ii2]; */
		}
	}

	/* 後退代入 */
	for(ii1=n-1; ii1>=0; ii1--){
		for(ii2=ii1+1; ii2<n; ii2++){
			SUB_MUL_COMPLEX(B[ii1], A[ii1][ii2], B[ii2]);
			/* B[ii1] -= A[ii1][ii2]*B[ii2]; */
		}
		B[ii1] = div_complex(B[ii1], A[ii1][ii1]);
		/* B[ii1] /= keep_away_zero(A[ii1][ii1]); */
	}

	return(NORMAL_RC);
}


/* a[n], b[n] の内積 */
COMPLEX inner_product_scalar_complex(int n, const double a[],
                                            const COMPLEX b[])
{
	int ii1;
	long double re, im;
	COMPLEX ret;

	re = 0.0;
	im = 0.0;
	for(ii1=0; ii1<n; ii1++){
		re += (long double)(a[ii1]) * (long double)(b[ii1].re);
		im += (long double)(a[ii1]) * (long double)(b[ii1].im);
	}

	ret.re = (double)re;
	ret.im = (double)im;

	return(ret);
}


/* a1[n], b[n] の内積, a2[n], b[n] の内積 */
RC inner_product_scalar2_complex(int n, const double a1[], const double a2[],
                                 const COMPLEX b[], COMPLEX *a1b, COMPLEX *a2b)
{
	int ii1;
	long double re1, im1;
	long double re2, im2;

	re1 = 0.0;
	im1 = 0.0;
	re2 = 0.0;
	im2 = 0.0;
	for(ii1=0; ii1<n; ii1++){
		re1 += (long double)(a1[ii1]) * (long double)(b[ii1].re);
		im1 += (long double)(a1[ii1]) * (long double)(b[ii1].im);
		re2 += (long double)(a2[ii1]) * (long double)(b[ii1].re);
		im2 += (long double)(a2[ii1]) * (long double)(b[ii1].im);
	}

	a1b->re = (double)re1;
	a1b->im = (double)im1;
	a2b->re = (double)re2;
	a2b->im = (double)im2;

	return(NORMAL_RC);
}


/* a1[n], b[n] の内積 => a1b */
/* ただし，添え字 nz_begin 〜 nz_end までを計算 */
RC inner_product_scalar1r_complex(int n, const double a1[],
                                  const COMPLEX b[], int nz_begin,
                                  int nz_end, COMPLEX *a1b)
{
	int ii1;
	long double re1, im1;

	re1 = im1 = 0.0;
	for(ii1=nz_begin; ii1<=nz_end; ii1++){
		re1 += (long double)(a1[ii1]) * (long double)(b[ii1].re);
		im1 += (long double)(a1[ii1]) * (long double)(b[ii1].im);
	}

	a1b->re = (double)re1;
	a1b->im = (double)im1;

	return(NORMAL_RC);
}


/* a1[n], b[n] の内積 => a1b */
/* a2[n], b[n] の内積 => a2b */
/* ただし，添え字 nz_begin 〜 nz_end までを計算 */
RC inner_product_scalar2r_complex(int n, const double a1[], const double a2[],
                                  const COMPLEX b[], int nz_begin, int nz_end,
                                  COMPLEX *a1b, COMPLEX *a2b)
{
	int ii1;
	long double re1, im1;
	long double re2, im2;

	re1 = im1 = 0.0;
	re2 = im2 = 0.0;
	for(ii1=nz_begin; ii1<=nz_end; ii1++){
		re1 += (long double)(a1[ii1]) * (long double)(b[ii1].re);
		im1 += (long double)(a1[ii1]) * (long double)(b[ii1].im);
		re2 += (long double)(a2[ii1]) * (long double)(b[ii1].re);
		im2 += (long double)(a2[ii1]) * (long double)(b[ii1].im);
	}

	a1b->re = (double)re1;
	a1b->im = (double)im1;
	a2b->re = (double)re2;
	a2b->im = (double)im2;

	return(NORMAL_RC);
}


/* a1[n], b[n] の内積 => a1b */
/* a2[n], b[n] の内積 => a2b */
/* a3[n], b[n] の内積 => a3b */
/* a4[n], b[n] の内積 => a4b */
/* a5[n], b[n] の内積 => a5b */
/* a6[n], b[n] の内積 => a6b */
/* ただし，添え字 nz_begin 〜 nz_end までを計算 */
RC inner_product_scalar6r_complex(int n, const double a1[], const double a2[],
                                  const double a3[], const double a4[],
                                  const double a5[], const double a6[],
                                  const COMPLEX b[], int nz_begin, int nz_end,
                                  COMPLEX *a1b, COMPLEX *a2b, COMPLEX *a3b,
                                  COMPLEX *a4b, COMPLEX *a5b, COMPLEX *a6b)
{
	int ii1;
	long double re1, im1;
	long double re2, im2;
	long double re3, im3;
	long double re4, im4;
	long double re5, im5;
	long double re6, im6;

	re1 = im1 = 0.0;
	re2 = im2 = 0.0;
	re3 = im3 = 0.0;
	re4 = im4 = 0.0;
	re5 = im5 = 0.0;
	re6 = im6 = 0.0;
	for(ii1=nz_begin; ii1<=nz_end; ii1++){
		re1 += (long double)(a1[ii1]) * (long double)(b[ii1].re);
		im1 += (long double)(a1[ii1]) * (long double)(b[ii1].im);
		re2 += (long double)(a2[ii1]) * (long double)(b[ii1].re);
		im2 += (long double)(a2[ii1]) * (long double)(b[ii1].im);
		re3 += (long double)(a3[ii1]) * (long double)(b[ii1].re);
		im3 += (long double)(a3[ii1]) * (long double)(b[ii1].im);
		re4 += (long double)(a4[ii1]) * (long double)(b[ii1].re);
		im4 += (long double)(a4[ii1]) * (long double)(b[ii1].im);
		re5 += (long double)(a5[ii1]) * (long double)(b[ii1].re);
		im5 += (long double)(a5[ii1]) * (long double)(b[ii1].im);
		re6 += (long double)(a6[ii1]) * (long double)(b[ii1].re);
		im6 += (long double)(a6[ii1]) * (long double)(b[ii1].im);
	}

	a1b->re = (double)re1;
	a1b->im = (double)im1;
	a2b->re = (double)re2;
	a2b->im = (double)im2;
	a3b->re = (double)re3;
	a3b->im = (double)im3;
	a4b->re = (double)re4;
	a4b->im = (double)im4;
	a5b->re = (double)re5;
	a5b->im = (double)im5;
	a6b->re = (double)re6;
	a6b->im = (double)im6;

	return(NORMAL_RC);
}


/* matrix[3][3]^T と v の積 -> v */
void mul_matrix33t_cvect(double matrix[][3], COMPLEX v[3])
{
	VECT3D re;
	VECT3D im;

	re.x = matrix[0][0]*v[0].re + matrix[1][0]*v[1].re + matrix[2][0]*v[2].re;
	re.y = matrix[0][1]*v[0].re + matrix[1][1]*v[1].re + matrix[2][1]*v[2].re;
	re.z = matrix[0][2]*v[0].re + matrix[1][2]*v[1].re + matrix[2][2]*v[2].re;

	im.x = matrix[0][0]*v[0].im + matrix[1][0]*v[1].im + matrix[2][0]*v[2].im;
	im.y = matrix[0][1]*v[0].im + matrix[1][1]*v[1].im + matrix[2][1]*v[2].im;
	im.z = matrix[0][2]*v[0].im + matrix[1][2]*v[1].im + matrix[2][2]*v[2].im;

	v[0].re = re.x;
	v[1].re = re.y;
	v[2].re = re.z;

	v[0].im = im.x;
	v[1].im = im.y;
	v[2].im = im.z;
}


