/*********************************************************************
 * nonzero_cg.h
 *
 * Copyright (C) 2001 AzLib Developers Group
 *
 * Written by
 *  <Kenzen TAKEUCHI> <Tomoyuki TSUBATA>
 *
 * Refer to log documents for details.
 *********************************************************************/

/* $Id: nonzero_cg.h 1059 2016-04-25 09:52:20Z hayashi $ */

#ifndef NONZERO_CG_H
#define NONZERO_CG_H

#include "rc.h"


typedef struct {
	int col;
	double v[3][3];
} NONZERO_ELEM3;

typedef struct {
	int size;
	NONZERO_ELEM3 *elem;
	int syn_size;
	int *syn;
} NONZERO_ROW3;

typedef struct {
	int size;
	NONZERO_ROW3 *row;
} NONZERO_MATRIX3;

typedef struct {
	int size;
	int diagonal_index;
	int *col;
	double *v;
} NONZERO_ROW;

typedef struct {
	int size;
	NONZERO_ROW *row;
} NONZERO_MATRIX;


size_t count_size_nonzero_matrix3(NONZERO_MATRIX3 matrix);
RC mul_nonzero_matrix_XtAX3s(NONZERO_MATRIX3 A, int X_size, double **X,
                             double **XtAX);
RC nonzero_B_orthogonalize3s(NONZERO_MATRIX3 B, int vect_size, 
                             double **vects, double **Bvects);
int count_nonzero_matrix3(NONZERO_MATRIX3 matrix);
RC delete_zero_matrix3s(NONZERO_MATRIX3 matrix);
RC print_nonzero_matrix3(FILE *fp, NONZERO_MATRIX3 matrix);
RC free_nonzero_matrix3(NONZERO_MATRIX3 *matrix);
RC iccg_nonzero3s(NONZERO_MATRIX3 matrix, const double *b, double *x,
                                       int precondition, int verb_sw);
RC modify_nonzero3s(NONZERO_MATRIX3 matrix,
                    int index, int xyz, double v, double *vect);
RC modify_nonzero3s_n(NONZERO_MATRIX3 matrix,
                      int index, int xyz, double v, double **vect,
                      int vect_size);
RC modal_const_nonzero3s(NONZERO_MATRIX3 *matrix,
                         int mode_size, const double mode[]);
NONZERO_ELEM3 *nonzero_matrix3s_ptr(NONZERO_MATRIX3 matrix, int row, int col);
double nonzero_matrix3s_val(NONZERO_MATRIX3 matrix, int row, int col);
RC row_col_zero_nonzero3s(NONZERO_MATRIX3 matrix, int index, int xyz);
RC row_col_print_nonzero3s(FILE *fp, NONZERO_MATRIX3 matrix,
                           int index, int xyz);
RC row_col_add_nonzero3s(NONZERO_MATRIX3 matrix, int index1, int xyz1,
                         double weight, int index2, int xyz2);
RC row_col_mul_nonzero3s(NONZERO_MATRIX3 matrix, int index,
                         double s_matrix[][3]);
RC nonzero3s_scaling_pre(NONZERO_MATRIX3 matrix,
                         double *scale_fact, double *vect);
RC nonzero3s_scaling_post(int matrix_size,
                          const double *scale_fact, double *vect);
RC mul_nonzero_matrix_vect3s(NONZERO_MATRIX3 matrix,
                             const double *vect, double *ans);


void print_nonzero_matrix1a(FILE *fp, NONZERO_MATRIX matrix);
void mul_scalar_nonzero_matrix1a(const double fact, NONZERO_MATRIX *matrix);
void mul_scalar_nonzero_row(const double fact, NONZERO_ROW *row);
void mul_nonzero_matrix_vect1a(NONZERO_MATRIX matrix, const double *vect,
                               double *ans);
void mul_nonzero_matrixT_vect1a(NONZERO_MATRIX matrix, const double *vect,
                                double *ans);
RC wadd_nonzero_matrix1a(double w1, NONZERO_MATRIX A, double w2,
		NONZERO_MATRIX B, NONZERO_MATRIX *C);
RC allocate_nonzero_matrix1a(int total_dof, NONZERO_MATRIX *matrix);
double *nonzero_row_ptr(NONZERO_ROW *n_row, int row, int col);
double *nonzero_matrix1a_ptr(NONZERO_MATRIX *matrix, int row, int col);
RC free_nonzero_matrix1a(NONZERO_MATRIX *matrix);
RC copy_nonzero_matrix1a(NONZERO_MATRIX source, NONZERO_MATRIX *copy);
RC copy_nonzero_matrix1a_offset(int offset_row, int offset_col, 
		NONZERO_MATRIX source, NONZERO_MATRIX *copy);
RC transpose_nonzero_matrix1a(NONZERO_MATRIX source, NONZERO_MATRIX *copy);
RC delete_zero_matrix1a (NONZERO_MATRIX *matrix);
RC nonzero1a_scaling(NONZERO_MATRIX *matrix, double *vector);
RC convert_nonzero1a2double(int row, int col, NONZERO_MATRIX source,
		double **dest);
RC convert_double2nonzero1a(int row, int col, double **source,
		NONZERO_MATRIX *dest);
RC cg_nonzero1a_ver1(NONZERO_MATRIX matrix, const double *b, double *x);
RC cg_nonzero1a_ver2(NONZERO_MATRIX matrix, const double *vector, double *ans);


#endif /* NONZERO_CG_H */


