/*********************************************************************
 * complex_utl.h
 *
 * Copyright (C) 2016 AzLib Developers Group
 *
 * Written by
 *  <Taiki AOYAMA> <Kohei SHINTANI> <Yuri NAKAMURA> 
 *  <Takuya HAYASHI>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: complex_utl.h 1094 2016-12-05 07:00:41Z hayashi $ */

#ifndef COMPLEX_UTL_H
#define COMPLEX_UTL_H

#include <complex.h>
#include "rc.h"
#include "nonzero_cg.h"

/* 虚数単位を表すマクロ I の再定義 */
#undef I
#ifdef _Imaginary_I
#define COMPLEX_UNIT (_Imaginary_I)
#else
#define COMPLEX_UNIT (_Complex_I)
#endif

typedef double _Complex COMPLEX;
#if 0 /* old */
typedef struct {
	double re;
	double im;
} COMPLEX;
#endif

/* 複素数3次元ベクトル（並進のみ） */
typedef struct {
	COMPLEX x;
	COMPLEX y;
	COMPLEX z;
} CVECT3D;

/* 複素数3次元ベクトル（並進 + 回転） */
typedef struct {
	COMPLEX x;
	COMPLEX y;
	COMPLEX z;
	COMPLEX xy;
	COMPLEX yz;
	COMPLEX zx;
} CVECT3DR;


/* CRS複素行 */
typedef struct {
	int size;
	int diagonal_index;
	int *col;
	COMPLEX *v;
} CRS_CROW;

/* CRS複素行列 */
typedef struct {
	int size;
	CRS_CROW *row;
} CRS_CMATRIX;


/* complex_utl.c */
COMPLEX polar2right_complex(double r, double theta);
RC allocate1C(int size, COMPLEX **array);
RC allocate2C(int record, int column, COMPLEX ***array);
RC allocate3C(int record1, int record2, int record3, COMPLEX ****array);
RC free1C(int size, COMPLEX **array);
RC free2C(int record, int column, COMPLEX ***array);
RC free3C(int record1, int record2, int record3, COMPLEX ****array);
RC copy_vect_complex(int size, const COMPLEX *vect1, COMPLEX *vect2);
RC init_vect_complex(int size, COMPLEX *vect);
RC init_matrix_complex(int row, int col, COMPLEX **matrix);
void mul_matrix_vect_complex (int row, int col, COMPLEX **matrix, 
		COMPLEX *vect, COMPLEX *ret);
RC print_vect_complex(FILE *fp, int size, const COMPLEX *vect);
RC print_matrix_complex(FILE *fp, int row, int col, const COMPLEX **matrix);
double norm_vect_complex(int size, COMPLEX *vect);
COMPLEX hermite_product(int size, COMPLEX *vect1, COMPLEX *vect2);

/* cvect3d.c */
RC conj_cvect3d(CVECT3D *vect);
double abs_cvect3d(CVECT3D vect);
double sq_abs_cvect3d(CVECT3D vect);
RC wadd_cvect3d(COMPLEX w1, CVECT3D vect1, COMPLEX w2, CVECT3D vect2, 
		CVECT3D *ret);
RC init_cvect3d(CVECT3D *vect);
RC init_cvect3dr(CVECT3DR *vect);
RC print_cvect3d(FILE *fp, CVECT3D vect);
RC print_cvect3dr(FILE *fp, CVECT3DR vect);
COMPLEX inner_product_cvect3d(CVECT3D vect1, CVECT3D vect2);
COMPLEX inner_product_cvect3dr(CVECT3DR vect1, CVECT3DR vect2);
RC tensor_product_cvect3d(CVECT3D vect1, CVECT3D vect2, COMPLEX **tensor);

/* crs_complex_utl.c */
RC mul_crs_cmatrix_cvect(CRS_CMATRIX matrix, COMPLEX *vect, 
        COMPLEX *ret);
RC mul_add_crs_cmatrix_cvect(CRS_CMATRIX matrix, COMPLEX *vect, 
        COMPLEX *ret);
RC copy_crs_cmatrix(CRS_CMATRIX source, CRS_CMATRIX *dest);
RC print_crs_cmatrix(FILE *fp, CRS_CMATRIX matrix);
RC print_crs_crow(FILE *fp, CRS_CROW row);
RC allocate_crs_cmatrix(int total_dof, CRS_CMATRIX *matrix);
RC free_crs_cmatrix(CRS_CMATRIX *matrix);
RC init_crs_cmatrix(CRS_CMATRIX *matrix);
COMPLEX *crs_crow_ptr(CRS_CROW *n_row, int row, int col);
COMPLEX *crs_cmatrix_ptr(CRS_CMATRIX *matrix, int row, int col);
RC wadd_crs_dmatrix2cmatrix(int part, double fact, NONZERO_MATRIX src, 
		CRS_CMATRIX *dst);
RC wadd_crs_drow2crow(int part, double fact, int row, NONZERO_ROW src,
		CRS_CROW *dst);
RC delete_zero_crs_cmatrix(CRS_CMATRIX *matrix);

#endif /* COMPLEX_UTL_H */

