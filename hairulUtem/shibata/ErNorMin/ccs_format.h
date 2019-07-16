#ifndef CRS_FORMAT_H
#define CRS_FORMAT_H

#include <stdio.h>
#include <stdlib.h>
//#include <float.h>
//#include <string.h>
#include "base.h"
#include "fem.h"
#include "mathematics.h"
  

typedef struct {
	int size;              /*   n  in PARDISO */
	int array_size;
	int alloc_array_size;
	int *row_offset;       /* ia[] in PARDISO */
	int *col_index;        /* ja[] in PARDISO */ 
	double *array;         /*  a[] in PARDISO */
} CRS;


typedef struct {
	int size;              /* n  in UMFPACK */
	int array_size;        /* nz in UMFPACK */
	int alloc_array_size;
	int *col_offset;       /* Ap[] in UMFPACK */
	int *row_index;        /* Ai[] in UMFPACK */ 
	int *col_index;        /* Tj[] in UMFPACK */ 
	double *array;         /* Ax[] in UMFPACK */
} CCS;

typedef struct {
	int size;
	int diagonal_index;
	int *row;
	double *v;
} NONZERO_COL;

typedef struct {
	int size;
	NONZERO_COL *col;
} NONZERO_MATRIX_CCS;


RC allocate_nonzero_ccs(int n, CCS *matrix);
RC clean_nonzero_ccs (CCS *matrix);
RC free_nonzero_ccs(CCS *matrix);
double *nonzero_ccs_ptr(CCS *matrix, int row, int col);
int search_nonzero_ccs_array_index(CCS matrix, int row, int col);
int search_nonzero_crs_array_index(CRS matrix, int row, int col);
RC modify_nonzero1a_crs(int dim, NODE_ARRAY node, BC_ARRAY rest,
                        NONZERO_MATRIX *matrix, double f_vect[]);
RC modify_nonzero1a_crs_dof(int index_xyz, double rest_value,
                                   NONZERO_MATRIX *matrix, double f_vect[]);
RC modify_ccs_format (int dim, NODE_ARRAY node, BC_ARRAY rest, CCS *ccs,
                      double f_vect[]);

/* Azlib like */
double *nonzero_matrix1a_ccs_ptr(NONZERO_MATRIX_CCS *matrix, int row, int col);
RC allocate_nonzero_matrix1a_ccs(int total_dof, NONZERO_MATRIX_CCS *matrix);
RC free_nonzero_matrix1a_ccs(NONZERO_MATRIX_CCS *matrix);
RC copy_nonzero_matrix1a_ccs2ccs_format(NONZERO_MATRIX_CCS *matrix, CCS *ccs);
RC copy_nonzero_matrix1a_crs2crs_format (NONZERO_MATRIX *matrix, CRS *crs);
RC free_nonzero_crs (CRS *matrix);

#endif
